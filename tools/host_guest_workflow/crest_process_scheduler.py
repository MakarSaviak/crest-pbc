#!/usr/bin/env python3
"""Deterministic process-isolated CREST ensemble scheduler.

Each frame is optimized by a fresh, single-threaded CREST process in a private
directory below --scratch-root.  Only compact validated results are copied to
--output.  The output directory is resumable: --resume skips only frames whose
record, input digest, output digest, and basic validation still agree.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import queue
import re
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

THREAD_ENV = {
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "BLIS_NUM_THREADS": "1",
    "OMP_DYNAMIC": "FALSE",
    "OMP_MAX_ACTIVE_LEVELS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
}
CALLS_RE = re.compile(r"Total number of energy\+grad calls:\s*(\d+)", re.I)
SUCCESS_RE = re.compile(r"(\d+)\s+of\s+(\d+)\s+structures successfully optimized", re.I)
FLOAT_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
ENERGY_RES = (
    re.compile(rf"(?:final\s+)?(?:single\s+point\s+)?energy\s*(?:=|:)?\s*({FLOAT_RE})", re.I),
    re.compile(rf"total\s+energy\s*(?:=|:)?\s*({FLOAT_RE})", re.I),
)
GRADIENT_RES = (
    re.compile(rf"(?:final\s+)?(?:rms|norm)\s+(?:of\s+)?(?:cartesian\s+)?gradient\s*(?:=|:)?\s*({FLOAT_RE})", re.I),
    re.compile(rf"(?:cartesian\s+)?gradient\s+(?:rms|norm)\s*(?:=|:)?\s*({FLOAT_RE})", re.I),
)
GUEST_GRADIENT_RES = (
    re.compile(rf"guest\s+(?:rms|norm)\s+(?:of\s+)?gradient\s*(?:=|:)?\s*({FLOAT_RE})", re.I),
    re.compile(rf"guest\s+gradient\s+(?:rms|norm)\s*(?:=|:)?\s*({FLOAT_RE})", re.I),
)


@dataclass
class FrameResult:
    frame: int
    worker: int | None
    attempt: int
    input_sha256: str
    start_s: float | None
    end_s: float | None
    wall_s: float | None
    calls: int | None
    converged: bool
    returncode: int | None
    final_energy_eh: float | None
    final_gradient_norm: float | None
    final_guest_gradient_norm: float | None
    gradient_source: str | None
    output_sha256: str | None
    host_max_displacement_a: float | None
    workspace_peak_bytes: int | None
    peak_rss_kib: int | None
    workspace_bytes_before_cleanup: int | None
    workspace_bytes_after_cleanup: int | None
    log_path: str | None
    output_path: str | None
    valid: bool
    error: str | None
    reused_from_resume: bool = False
    reference_calls_match: bool | None = None
    reference_energy_error_eh: float | None = None
    reference_guest_rmsd_a: float | None = None
    reference_output_sha256_match: bool | None = None


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(value)
    os.replace(temporary, path)


def atomic_write_bytes(path: Path, value: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_bytes(value)
    os.replace(temporary, path)


def split_xyz(path: Path) -> list[str]:
    lines = path.read_text().splitlines()
    frames: list[str] = []
    cursor = 0
    while cursor < len(lines):
        if not lines[cursor].strip():
            cursor += 1
            continue
        try:
            natoms = int(lines[cursor].strip())
        except ValueError as exc:
            raise ValueError(f"{path}: invalid XYZ atom count at line {cursor + 1}") from exc
        stop = cursor + natoms + 2
        if natoms <= 0 or stop > len(lines):
            raise ValueError(f"{path}: truncated XYZ frame at line {cursor + 1}")
        frames.append("\n".join(lines[cursor:stop]) + "\n")
        cursor = stop
    if not frames:
        raise ValueError(f"{path}: contains no XYZ frames")
    return frames


def xyz_coordinates(frame: str) -> list[tuple[float, float, float]]:
    lines = frame.splitlines()
    natoms = int(lines[0])
    if len(lines) != natoms + 2:
        raise ValueError("invalid XYZ frame")
    values: list[tuple[float, float, float]] = []
    for line in lines[2:]:
        fields = line.split()
        if len(fields) < 4:
            raise ValueError(f"invalid XYZ coordinate line: {line!r}")
        values.append((float(fields[1]), float(fields[2]), float(fields[3])))
    return values


def comment_energy(frame: str) -> float | None:
    comment = frame.splitlines()[1].strip()
    first = comment.split(maxsplit=1)[0] if comment else ""
    try:
        return float(first.replace("D", "E").replace("d", "e"))
    except ValueError:
        for pattern in ENERGY_RES:
            found = pattern.search(comment)
            if found:
                return float(found.group(1).replace("D", "E").replace("d", "e"))
    return None


def parse_log_energy(text: str) -> float | None:
    matches: list[float] = []
    for pattern in ENERGY_RES:
        matches.extend(float(m.group(1).replace("D", "E").replace("d", "e")) for m in pattern.finditer(text))
    return matches[-1] if matches else None


def parse_log_gradient_norm(text: str) -> float | None:
    matches: list[float] = []
    for pattern in GRADIENT_RES:
        matches.extend(float(m.group(1).replace("D", "E").replace("d", "e")) for m in pattern.finditer(text))
    return matches[-1] if matches else None


def parse_log_guest_gradient_norm(text: str) -> float | None:
    matches: list[float] = []
    for pattern in GUEST_GRADIENT_RES:
        matches.extend(float(m.group(1).replace("D", "E").replace("d", "e")) for m in pattern.finditer(text))
    return matches[-1] if matches else None


def parse_calls(text: str) -> int | None:
    match = CALLS_RE.search(text)
    return int(match.group(1)) if match else None


def parse_converged(text: str, returncode: int, nframes: int) -> bool:
    match = SUCCESS_RE.search(text)
    return returncode == 0 and match is not None and int(match.group(1)) == nframes


def rmsd_direct(left: list[tuple[float, float, float]], right: list[tuple[float, float, float]]) -> float:
    if len(left) != len(right):
        raise ValueError("coordinate count differs")
    if not left:
        return 0.0
    return math.sqrt(
        sum((a - b) ** 2 for xyz1, xyz2 in zip(left, right) for a, b in zip(xyz1, xyz2))
        / len(left)
    )


def host_max_displacement(
    initial: list[tuple[float, float, float]], final: list[tuple[float, float, float]], frozen_atoms: int
) -> float:
    if len(initial) != len(final):
        raise ValueError("coordinate count differs")
    if not 0 <= frozen_atoms <= len(initial):
        raise ValueError("invalid frozen atom count")
    return max(
        (math.sqrt(sum((a - b) ** 2 for a, b in zip(initial[i], final[i]))) for i in range(frozen_atoms)),
        default=0.0,
    )


def dir_size(path: Path) -> int:
    total = 0
    if not path.exists():
        return total
    for candidate in path.rglob("*"):
        try:
            if candidate.is_file():
                total += candidate.stat().st_size
        except FileNotFoundError:
            pass
    return total


class WorkspaceMonitor:
    def __init__(self, path: Path, interval_s: float) -> None:
        self.path = path
        self.interval_s = interval_s
        self.peak_bytes = 0
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._sample, daemon=True)

    def _sample(self) -> None:
        while not self._stop.is_set():
            self.peak_bytes = max(self.peak_bytes, dir_size(self.path))
            self._stop.wait(self.interval_s)

    def start(self) -> None:
        self._thread.start()

    def stop(self) -> int:
        self._stop.set()
        self._thread.join()
        self.peak_bytes = max(self.peak_bytes, dir_size(self.path))
        return self.peak_bytes


def clean_workspace(path: Path, keep_workspaces: bool) -> int:
    before = dir_size(path)
    if not keep_workspaces:
        shutil.rmtree(path, ignore_errors=True)
    return before if keep_workspaces else dir_size(path)


def command_for_frame(args: argparse.Namespace, workspace: Path, resource_file: Path) -> list[str]:
    command = [
        args.time_exe,
        "-f",
        "%M",
        "-o",
        str(resource_file),
        args.exe,
        "--mdopt",
        "input.xyz",
        "--gfnff",
        "--cinp",
        "constraints.inp",
        "-T",
        "1",
        "-opt",
        args.opt_level,
    ]
    command.extend(args.crest_arg)
    return command


def read_peak_rss_kib(path: Path) -> int | None:
    try:
        return int(path.read_text().strip().splitlines()[-1])
    except (FileNotFoundError, IndexError, ValueError):
        return None


def result_path(output: Path, frame: int) -> Path:
    return output / "results" / f"f{frame:05d}.json"


def frame_path(output: Path, frame: int) -> Path:
    return output / "results" / f"f{frame:05d}.xyz"


def log_path(output: Path, frame: int) -> Path:
    return output / "logs" / f"f{frame:05d}.out"


def load_existing_result(output: Path, frame: int, input_sha: str) -> FrameResult | None:
    record = result_path(output, frame)
    final_xyz = frame_path(output, frame)
    if not record.exists() or not final_xyz.exists():
        return None
    try:
        raw = json.loads(record.read_text())
        if not raw.get("valid", False) or not raw.get("converged", False):
            return None
        if raw.get("input_sha256") != input_sha:
            return None
        if raw.get("output_sha256") != sha256_file(final_xyz):
            return None
        allowed = set(FrameResult.__dataclass_fields__)
        raw = {key: value for key, value in raw.items() if key in allowed}
        raw["reused_from_resume"] = True
        return FrameResult(**raw)
    except (OSError, ValueError, TypeError):
        return None


def load_costs(path: Path | None) -> dict[int, float]:
    if path is None:
        return {}
    if path.suffix.lower() == ".json":
        raw = json.loads(path.read_text())
        if isinstance(raw, dict) and "frames" in raw:
            raw = raw["frames"]
        if isinstance(raw, list):
            values = {int(item["frame"]): float(item.get("predicted_cost", item.get("calls", item.get("wall_s")))) for item in raw}
        else:
            values = {int(key): float(value) for key, value in raw.items()}
        return values
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    values = {}
    for row in rows:
        score = row.get("predicted_cost") or row.get("predicted_cost_s") or row.get("calls") or row.get("wall_s")
        if score is None:
            raise ValueError(f"{path}: expected predicted_cost, calls, or wall_s column")
        values[int(row["frame"])] = float(score)
    return values


def queue_order(frame_ids: list[int], args: argparse.Namespace, costs: dict[int, float]) -> list[int]:
    if args.order == "fifo":
        return frame_ids
    if args.order == "predicted-longest-first":
        return sorted(frame_ids, key=lambda frame: (-costs.get(frame, float("-inf")), frame))
    raise ValueError(args.order)


def run_one_attempt(
    args: argparse.Namespace,
    frame_id: int,
    frame: str,
    input_sha: str,
    worker_id: int,
    attempt: int,
    scratch_run: Path,
    global_start: float,
) -> FrameResult:
    workspace = scratch_run / f"f{frame_id:05d}" / f"attempt-{attempt:02d}"
    workspace.mkdir(parents=True, exist_ok=False)
    (workspace / "input.xyz").write_text(frame)
    shutil.copy2(args.constraints, workspace / "constraints.inp")
    env = os.environ.copy()
    env.update(THREAD_ENV)
    resource = workspace / "resource.maxrss_kib"
    run_log = workspace / "crest.out"
    monitor = WorkspaceMonitor(workspace, args.disk_sample_interval)
    start = time.perf_counter()
    completed: subprocess.CompletedProcess[str] | None = None
    monitor.start()
    try:
        with run_log.open("w") as handle:
            completed = subprocess.run(
                command_for_frame(args, workspace, resource),
                cwd=workspace,
                env=env,
                stdout=handle,
                stderr=subprocess.STDOUT,
                check=False,
            )
        end = time.perf_counter()
        text = run_log.read_text(errors="replace")
        output = workspace / "crest_ensemble.xyz"
        converged = parse_converged(text, completed.returncode, 1)
        if not output.exists():
            raise RuntimeError("CREST did not write crest_ensemble.xyz")
        output_text = output.read_text()
        output_coords = xyz_coordinates(output_text)
        initial_coords = xyz_coordinates(frame)
        if len(output_coords) != len(initial_coords):
            raise RuntimeError(f"CREST wrote {len(output_coords)} atoms; expected {len(initial_coords)}")
        output_energy = comment_energy(output_text)
        energy = output_energy if output_energy is not None else parse_log_energy(text)
        gradient = parse_log_gradient_norm(text)
        guest_gradient = parse_log_guest_gradient_norm(text)
        host_displacement = host_max_displacement(initial_coords, output_coords, args.frozen_atoms)
        before_cleanup = dir_size(workspace)
        peak = monitor.stop()
        peak_rss = read_peak_rss_kib(resource)
        final_xyz = frame_path(args.output, frame_id)
        final_log = log_path(args.output, frame_id)
        atomic_write_bytes(final_xyz, output.read_bytes())
        atomic_write_text(final_log, text)
        valid = converged and host_displacement == 0.0
        result = FrameResult(
            frame=frame_id,
            worker=worker_id,
            attempt=attempt,
            input_sha256=input_sha,
            start_s=start - global_start,
            end_s=end - global_start,
            wall_s=end - start,
            calls=parse_calls(text),
            converged=converged,
            returncode=completed.returncode,
            final_energy_eh=energy,
            final_gradient_norm=gradient,
            final_guest_gradient_norm=guest_gradient,
            gradient_source=(
                "guest-specific CREST log field"
                if guest_gradient is not None
                else ("CREST log only reports a non-guest gradient norm" if gradient is not None else None)
            ),
            output_sha256=sha256_file(final_xyz),
            host_max_displacement_a=host_displacement,
            workspace_peak_bytes=peak,
            peak_rss_kib=peak_rss,
            workspace_bytes_before_cleanup=before_cleanup,
            workspace_bytes_after_cleanup=clean_workspace(workspace, args.keep_workspaces),
            log_path=str(final_log.relative_to(args.output)),
            output_path=str(final_xyz.relative_to(args.output)),
            valid=valid,
            error=None if valid else "CREST did not report one converged structure or host coordinates changed",
        )
    except Exception as exc:
        end = time.perf_counter()
        peak = monitor.stop()
        text = run_log.read_text(errors="replace") if run_log.exists() else ""
        before_cleanup = dir_size(workspace)
        final_log = log_path(args.output, frame_id)
        if text:
            atomic_write_text(final_log, text)
        result = FrameResult(
            frame=frame_id,
            worker=worker_id,
            attempt=attempt,
            input_sha256=input_sha,
            start_s=start - global_start,
            end_s=end - global_start,
            wall_s=end - start,
            calls=parse_calls(text),
            converged=False,
            returncode=completed.returncode if completed is not None else None,
            final_energy_eh=parse_log_energy(text),
            final_gradient_norm=parse_log_gradient_norm(text),
            final_guest_gradient_norm=parse_log_guest_gradient_norm(text),
            gradient_source=None,
            output_sha256=None,
            host_max_displacement_a=None,
            workspace_peak_bytes=peak,
            peak_rss_kib=read_peak_rss_kib(resource),
            workspace_bytes_before_cleanup=before_cleanup,
            workspace_bytes_after_cleanup=clean_workspace(workspace, args.keep_workspaces),
            log_path=str(final_log.relative_to(args.output)) if final_log.exists() else None,
            output_path=None,
            valid=False,
            error=str(exc),
        )
    atomic_write_text(result_path(args.output, frame_id), json.dumps(asdict(result), indent=2, sort_keys=True) + "\n")
    return result


def run_frame(
    args: argparse.Namespace,
    frame_id: int,
    frame: str,
    worker_id: int,
    scratch_run: Path,
    global_start: float,
) -> FrameResult:
    input_sha = sha256_bytes(frame.encode())
    existing = load_existing_result(args.output, frame_id, input_sha) if args.resume else None
    if existing is not None:
        return existing
    frame_path(args.output, frame_id).unlink(missing_ok=True)
    result_path(args.output, frame_id).unlink(missing_ok=True)
    last: FrameResult | None = None
    for attempt in range(1, args.retries + 2):
        last = run_one_attempt(args, frame_id, frame, input_sha, worker_id, attempt, scratch_run, global_start)
        if last.valid:
            return last
    assert last is not None
    return last


def reference_row(reference: Path, frame: int) -> dict[str, str] | None:
    table = reference / "frames.csv"
    if not table.exists():
        return None
    with table.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if int(row["frame"]) == frame:
                return row
    return None


def validate_against_reference(result: FrameResult, args: argparse.Namespace, initial_frame: str) -> FrameResult:
    if args.reference is None or not result.valid:
        return result
    reference_xyz = frame_path(args.reference, result.frame)
    if not reference_xyz.exists():
        result.valid = False
        result.error = f"reference output missing: {reference_xyz}"
        return result
    ours = frame_path(args.output, result.frame).read_text()
    ref = reference_xyz.read_text()
    ours_coords = xyz_coordinates(ours)
    ref_coords = xyz_coordinates(ref)
    result.reference_guest_rmsd_a = rmsd_direct(ours_coords[args.frozen_atoms :], ref_coords[args.frozen_atoms :])
    ours_energy = comment_energy(ours)
    ref_energy = comment_energy(ref)
    if ours_energy is not None and ref_energy is not None:
        result.reference_energy_error_eh = abs(ours_energy - ref_energy)
    result.reference_output_sha256_match = sha256_bytes(ours.encode()) == sha256_bytes(ref.encode())
    ref_row = reference_row(args.reference, result.frame)
    if ref_row and ref_row.get("calls"):
        result.reference_calls_match = result.calls == int(ref_row["calls"])
    failed_checks = []
    if result.reference_calls_match is False:
        failed_checks.append("energy-gradient call count differs from reference")
    if result.reference_energy_error_eh is not None and result.reference_energy_error_eh > args.energy_tolerance_eh:
        failed_checks.append("final energy differs from reference")
    if result.reference_guest_rmsd_a > args.guest_rmsd_tolerance_a:
        failed_checks.append("guest RMSD differs from reference")
    if args.require_byte_identical and not result.reference_output_sha256_match:
        failed_checks.append("output SHA-256 differs from reference")
    if failed_checks:
        result.valid = False
        result.error = "; ".join(failed_checks)
    atomic_write_text(result_path(args.output, result.frame), json.dumps(asdict(result), indent=2, sort_keys=True) + "\n")
    return result


def worker_assignments(frame_ids: list[int], workers: int, schedule: str) -> list[list[int]]:
    assignments = [[] for _ in range(workers)]
    if schedule == "static-contiguous":
        width = math.ceil(len(frame_ids) / workers)
        for worker in range(workers):
            assignments[worker] = frame_ids[worker * width : (worker + 1) * width]
    elif schedule == "static-roundrobin":
        for index, frame_id in enumerate(frame_ids):
            assignments[index % workers].append(frame_id)
    else:
        raise ValueError(schedule)
    return assignments


def run_scheduler(args: argparse.Namespace) -> int:
    args.input = args.input.resolve()
    args.constraints = args.constraints.resolve()
    args.output = args.output.resolve()
    args.scratch_root = args.scratch_root.resolve()
    args.exe = str(Path(args.exe).resolve())
    args.time_exe = str(Path(args.time_exe).resolve())
    if args.reference is not None:
        args.reference = args.reference.resolve()
    if args.predicted_costs is not None:
        args.predicted_costs = args.predicted_costs.resolve()
    frames = split_xyz(args.input)
    if not 0 <= args.frozen_atoms < len(xyz_coordinates(frames[0])):
        raise ValueError("--frozen-atoms must be non-negative and smaller than the atom count")
    if args.workers < 1:
        raise ValueError("--workers must be positive")
    if args.workers > len(frames):
        print(f"warning: capping workers from {args.workers} to {len(frames)} frames", file=sys.stderr)
        args.workers = len(frames)
    prepare_output(args)
    costs = load_costs(args.predicted_costs)
    input_digests = [sha256_bytes(frame.encode()) for frame in frames]
    scratch_run = args.scratch_root / args.run_id
    if scratch_run.exists() and not args.resume:
        raise FileExistsError(f"scratch run directory already exists: {scratch_run}; choose --run-id or use --resume")
    scratch_run.mkdir(parents=True, exist_ok=True)
    manifest = {
        "schema": 3,
        "input": str(args.input.resolve()),
        "input_sha256": sha256_file(args.input),
        "frame_input_sha256": input_digests,
        "constraints_sha256": sha256_file(args.constraints),
        "crest_exe": str(Path(args.exe).resolve()),
        "crest_exe_sha256": sha256_file(Path(args.exe)),
        "crest_args": args.crest_arg,
        "opt_level": args.opt_level,
        "frozen_atoms": args.frozen_atoms,
        "thread_environment": THREAD_ENV,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    check_or_write_manifest(args, manifest)
    all_ids = list(range(len(frames)))
    ordered = queue_order(all_ids, args, costs)
    existing: dict[int, FrameResult] = {}
    if args.resume:
        for frame_id, input_sha in enumerate(input_digests):
            loaded = load_existing_result(args.output, frame_id, input_sha)
            if loaded is not None:
                loaded = validate_against_reference(loaded, args, frames[frame_id])
                if loaded.valid:
                    existing[frame_id] = loaded
    pending = [frame_id for frame_id in ordered if frame_id not in existing]
    global_start = time.perf_counter()
    results = list(existing.values())
    result_lock = threading.Lock()

    def complete(result: FrameResult) -> None:
        checked = validate_against_reference(result, args, frames[result.frame])
        with result_lock:
            results.append(checked)
        label = "resume" if checked.reused_from_resume else f"worker={checked.worker}"
        print(
            f"frame={checked.frame} {label} attempt={checked.attempt} "
            f"calls={checked.calls} wall_s={checked.wall_s} valid={checked.valid}",
            flush=True,
        )

    if args.schedule == "dynamic":
        task_queue: queue.Queue[int] = queue.Queue()
        for frame_id in pending:
            task_queue.put(frame_id)

        def dynamic_worker(worker_id: int) -> None:
            while True:
                try:
                    frame_id = task_queue.get_nowait()
                except queue.Empty:
                    return
                try:
                    complete(run_frame(args, frame_id, frames[frame_id], worker_id, scratch_run, global_start))
                finally:
                    task_queue.task_done()

        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(dynamic_worker, worker_id) for worker_id in range(args.workers)]
            for future in as_completed(futures):
                future.result()
    else:
        assignments = worker_assignments(pending, args.workers, args.schedule)

        def static_worker(worker_id: int, frame_ids: list[int]) -> None:
            for frame_id in frame_ids:
                complete(run_frame(args, frame_id, frames[frame_id], worker_id, scratch_run, global_start))

        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(static_worker, worker, assigned) for worker, assigned in enumerate(assignments)]
            for future in as_completed(futures):
                future.result()

    global_end = time.perf_counter()
    results.sort(key=lambda record: record.frame)
    write_outputs(args, frames, results, ordered, costs, global_end - global_start)
    if not args.keep_workspaces:
        shutil.rmtree(scratch_run, ignore_errors=True)
    failures = [record for record in results if not record.valid]
    return 1 if failures or len(results) != len(frames) else 0


def prepare_output(args: argparse.Namespace) -> None:
    if args.output.exists() and any(args.output.iterdir()) and not args.resume and not args.overwrite:
        raise FileExistsError(f"{args.output} is non-empty; use --resume or --overwrite")
    if args.overwrite:
        if args.output.exists():
            shutil.rmtree(args.output)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "results").mkdir(exist_ok=True)
    (args.output / "logs").mkdir(exist_ok=True)


def check_or_write_manifest(args: argparse.Namespace, manifest: dict[str, Any]) -> None:
    path = args.output / "manifest.json"
    stable = (
        "schema",
        "input_sha256",
        "frame_input_sha256",
        "constraints_sha256",
        "crest_exe_sha256",
        "crest_args",
        "opt_level",
        "frozen_atoms",
        "thread_environment",
    )
    if path.exists() and args.resume:
        previous = json.loads(path.read_text())
        changed = [key for key in stable if previous.get(key) != manifest.get(key)]
        if changed:
            raise ValueError(f"resume manifest differs for: {', '.join(changed)}")
        return
    atomic_write_text(path, json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def write_outputs(
    args: argparse.Namespace,
    frames: list[str],
    results: list[FrameResult],
    ordered: list[int],
    costs: dict[int, float],
    makespan: float,
) -> None:
    records = [asdict(result) for result in results]
    with (args.output / "frames.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(FrameResult.__dataclass_fields__))
        writer.writeheader()
        writer.writerows(records)
    atomic_write_text(args.output / "frames.json", json.dumps(records, indent=2, sort_keys=True) + "\n")
    successful = [result for result in results if result.valid]
    if len(successful) == len(frames):
        ensemble = b"".join(frame_path(args.output, frame).read_bytes() for frame in range(len(frames)))
        atomic_write_bytes(args.output / "crest_ensemble.xyz", ensemble)
    else:
        ensemble = b""
    resumed_frames = sum(result.reused_from_resume for result in results)
    invocation_busy = sum(result.wall_s or 0.0 for result in results if not result.reused_from_resume)
    cumulative_frame_busy = sum(result.wall_s or 0.0 for result in results)
    worker_end: dict[int, float] = {worker: 0.0 for worker in range(args.workers)}
    worker_busy: dict[int, float] = {worker: 0.0 for worker in range(args.workers)}
    for result in results:
        if result.worker is not None and not result.reused_from_resume:
            worker_end[result.worker] = max(worker_end[result.worker], result.end_s or 0.0)
            worker_busy[result.worker] += result.wall_s or 0.0
    used_end = [value for value in worker_end.values() if value > 0.0]
    tail_idle_core_s = sum(max(used_end) - value for value in used_end) if used_end else 0.0
    # Frame timestamps from prior invocations have a different monotonic-time
    # origin.  A final resume segment must therefore never be labelled as the
    # end-to-end makespan.  Slurm accounting (or another orchestration ledger)
    # is required to combine separate allocations faithfully.
    timing_complete = resumed_frames == 0
    summary = {
        "schema": 3,
        "frames": len(frames),
        "workers": args.workers,
        "schedule": args.schedule,
        "order": args.order,
        "queue_order": ordered,
        "predicted_cost": {str(key): costs.get(key) for key in ordered},
        "retries": args.retries,
        "resumed_frames": resumed_frames,
        "completed_frames": len(successful),
        "failed_frames": len(frames) - len(successful),
        "failure_attempts": sum(max(0, result.attempt - 1) for result in results),
        "timing_status": "single_invocation" if timing_complete else "resume_segment_only",
        "makespan_s": makespan if timing_complete else None,
        "throughput_frames_s": len(successful) / makespan if timing_complete and makespan else None,
        "aggregate_busy_s": invocation_busy,
        "cumulative_frame_busy_s": cumulative_frame_busy,
        "invocation_makespan_s": makespan,
        "requested_core_utilization": invocation_busy / (args.workers * makespan) if timing_complete and makespan else None,
        "worker_finish_spread_s": max(used_end) - min(used_end) if used_end else 0.0,
        "tail_idle_core_s": tail_idle_core_s,
        "worker_busy_s": worker_busy,
        "worker_end_s": worker_end,
        "peak_workspace_bytes": max((result.workspace_peak_bytes or 0 for result in results), default=0),
        "aggregate_workspace_peak_bytes": sum(result.workspace_peak_bytes or 0 for result in results),
        "peak_rss_kib": max((result.peak_rss_kib or 0 for result in results), default=0),
        "retained_workspace_bytes": sum(result.workspace_bytes_after_cleanup or 0 for result in results),
        "total_calls": sum(result.calls or 0 for result in results),
        "ensemble_sha256": sha256_bytes(ensemble) if ensemble else None,
        "validation_passed": len(successful) == len(frames),
    }
    atomic_write_text(args.output / "summary.json", json.dumps(summary, indent=2, sort_keys=True) + "\n")


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--input", type=Path, required=True)
    value.add_argument("--constraints", type=Path, required=True)
    value.add_argument("--exe", required=True, help="CREST executable")
    value.add_argument("--output", type=Path, required=True)
    value.add_argument("--scratch-root", type=Path, required=True)
    value.add_argument("--run-id", default="crest-queue")
    value.add_argument("--workers", type=int, required=True)
    value.add_argument("--schedule", choices=("dynamic", "static-contiguous", "static-roundrobin"), default="dynamic")
    value.add_argument("--order", choices=("fifo", "predicted-longest-first"), default="fifo")
    value.add_argument("--predicted-costs", type=Path)
    value.add_argument("--retries", type=int, default=0)
    value.add_argument("--resume", action="store_true")
    value.add_argument("--overwrite", action="store_true")
    value.add_argument("--keep-workspaces", action="store_true")
    value.add_argument("--frozen-atoms", type=int, default=0, help="number of leading frozen atoms; set explicitly for frozen-host jobs")
    value.add_argument("--opt-level", default="vloose")
    value.add_argument("--crest-arg", action="append", default=[], help="extra argument passed to CREST; repeatable")
    value.add_argument("--time-exe", default="/usr/bin/time")
    value.add_argument("--disk-sample-interval", type=float, default=0.25)
    value.add_argument("--reference", type=Path, help="v3 one-core reference output directory")
    value.add_argument("--energy-tolerance-eh", type=float, default=1.0e-8)
    value.add_argument("--guest-rmsd-tolerance-a", type=float, default=1.0e-6)
    value.add_argument("--require-byte-identical", action="store_true")
    return value


if __name__ == "__main__":
    try:
        raise SystemExit(run_scheduler(parser().parse_args()))
    except Exception as failure:
        print(f"crest_queue_scheduler_v3: {failure}", file=sys.stderr)
        raise SystemExit(2)
