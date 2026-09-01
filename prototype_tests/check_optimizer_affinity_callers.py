#!/usr/bin/env python3
"""Fail closed unless every optimizer caller returns before continuation work."""

from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path


DIRECT_EXPECTED = Counter(
    {
        "crest_search_1": 1,
        "crest_ensemble_optimization": 1,
        "crest_refine": 1,
        "crest_search_mecp": 1,
        "crest_new_protonate": 3,
        "crest_new_deprotonate": 3,
        "crest_new_tautomerize": 3,
        "crest_multilevel_oloop": 1,
    }
)
REFINE_EXPECTED = Counter(
    {
        "crest_ensemble_optimization": 1,
        "crest_new_protonate": 1,
        "crest_new_deprotonate": 1,
        "crest_new_tautomerize": 1,
        "crest_multilevel_oloop": 1,
    }
)

SUBROUTINE_RE = re.compile(r"^\s*subroutine\s+([a-z0-9_]+)", re.IGNORECASE)
INTERFACE_RE = re.compile(r"^\s*(?:abstract\s+)?interface(?:\s|$)", re.IGNORECASE)
END_INTERFACE_RE = re.compile(r"^\s*end\s*interface(?:\s|$)", re.IGNORECASE)
DIRECT_RE = re.compile(r"^\s*call\s+crest_oloop\s*\(", re.IGNORECASE)
REFINE_RE = re.compile(r"^\s*call\s+crest_refine\s*\(", re.IGNORECASE)
INLINE_REFINE_RE = re.compile(
    r"^\s*if\s*\(.*\)\s*call\s+crest_refine\s*\(", re.IGNORECASE
)
GUARD_RE = re.compile(
    r"^\s*if\s*\(\s*env%iostatus_meta\s*/=\s*status_normal\s*\)\s*then\s*$",
    re.IGNORECASE,
)
BLOCK_IF_RE = re.compile(r"^\s*if\s*\(.*\)\s*then\s*$", re.IGNORECASE)
END_IF_RE = re.compile(r"^\s*end\s*if\s*$|^\s*endif\s*$", re.IGNORECASE)
ELSE_RE = re.compile(r"^\s*else\s*$", re.IGNORECASE)
RETURN_RE = re.compile(r"^\s*return\s*$", re.IGNORECASE)
GOTO_RE = re.compile(r"^\s*goto\s+([0-9]+)\s*$", re.IGNORECASE)
LABEL_CONTINUE_RE = re.compile(r"^\s*([0-9]+)\s+continue\b", re.IGNORECASE)
END_SUBROUTINE_RE = re.compile(r"^\s*end\s+subroutine\b", re.IGNORECASE)
SAFE_CLEANUP_RES = (
    re.compile(r"^\s*call\s+tim%stop\s*\(", re.IGNORECASE),
    re.compile(r"^\s*deallocate\s*\(", re.IGNORECASE),
    re.compile(r"^\s*env%[a-z0-9_%]+\s*=", re.IGNORECASE),
)


def significant(line: str) -> bool:
    stripped = line.strip()
    return bool(stripped) and not stripped.startswith("!")


def next_significant(lines: list[str], start: int) -> int:
    for index in range(start, len(lines)):
        if significant(lines[index]):
            return index
    raise ValueError("call has no following executable statement")


def statement_end(lines: list[str], start: int) -> int:
    """Return the last physical line of a continued Fortran statement."""
    index = start
    while lines[index].rstrip().endswith("&"):
        index = next_significant(lines, index + 1)
    return index


def guard_after_call(
    lines: list[str], call_index: int, kind: str, inline_refine: bool, location: str
) -> tuple[int, set[int]]:
    """Locate the shared typed-status guard after one logical optimizer call.

    Most callers are a call followed immediately by the guard.  The
    multilevel caller has two alternative, continued calls in an `if/else`
    dispatch and a single guard after `end if`; treat that as one semantic
    call only after proving that each branch consists solely of the call.
    """
    call_end = statement_end(lines, call_index)
    next_index = next_significant(lines, call_end + 1)
    if GUARD_RE.match(lines[next_index]):
        return next_index, set()
    if inline_refine and "Poststage refinement wall time:" in lines[next_index]:
        timing_end = statement_end(lines, next_index)
        guard_index = next_significant(lines, timing_end + 1)
        if GUARD_RE.match(lines[guard_index]):
            return guard_index, set()
        raise ValueError(
            f"{location}: timing observation is not followed by the typed-status guard"
        )
    if not ELSE_RE.match(lines[next_index]):
        raise ValueError(
            f"{location}: first executable statement after call is not the typed-status guard"
        )

    alternate_index = next_significant(lines, next_index + 1)
    alternate_kind = (
        "direct" if DIRECT_RE.match(lines[alternate_index]) else
        "refine" if REFINE_RE.match(lines[alternate_index]) else ""
    )
    if alternate_kind != kind:
        raise ValueError(
            f"{location}: alternate branch does not contain the same optimizer call"
        )
    alternate_end = statement_end(lines, alternate_index)
    end_if_index = next_significant(lines, alternate_end + 1)
    if not END_IF_RE.match(lines[end_if_index]):
        raise ValueError(
            f"{location}: alternate optimizer branch contains continuation work"
        )
    guard_index = next_significant(lines, end_if_index + 1)
    if not GUARD_RE.match(lines[guard_index]):
        raise ValueError(
            f"{location}: shared branch guard is not the typed-status guard"
        )
    return guard_index, {alternate_index}


def goto_reaches_terminal_cleanup(lines: list[str], start: int, label: str) -> bool:
    """Accept only a forward jump to cleanup that returns before another optimizer call."""
    target = -1
    for index in range(start + 1, len(lines)):
        match = LABEL_CONTINUE_RE.match(lines[index])
        if match and match.group(1) == label:
            target = index
            break
        if END_SUBROUTINE_RE.match(lines[index]):
            return False
    if target < 0:
        return False
    for index in range(target + 1, len(lines)):
        if DIRECT_RE.match(lines[index]) or REFINE_RE.match(lines[index]):
            return False
        if RETURN_RE.match(lines[index]):
            return True
        if END_SUBROUTINE_RE.match(lines[index]):
            return False
    return False


def validate_guard(lines: list[str], guard_index: int, location: str) -> None:
    if not GUARD_RE.match(lines[guard_index]):
        raise ValueError(
            f"{location}: first executable statement after call is not the typed-status guard"
        )

    depth = 0
    return_seen = False
    end_index = -1
    for index in range(guard_index, len(lines)):
        line = lines[index]
        if not significant(line):
            continue
        if BLOCK_IF_RE.match(line):
            depth += 1
            continue
        if END_IF_RE.match(line):
            depth -= 1
            if depth == 0:
                end_index = index
                break
            if depth < 0:
                raise ValueError(f"{location}: malformed failure guard")
            continue
        if index == guard_index:
            continue
        if RETURN_RE.match(line):
            if depth != 1:
                raise ValueError(f"{location}: return is not in the immediate failure branch")
            return_seen = True
            continue
        goto_match = GOTO_RE.match(line)
        if goto_match:
            if depth != 1:
                raise ValueError(f"{location}: goto is not in the immediate failure branch")
            if not goto_reaches_terminal_cleanup(lines, index, goto_match.group(1)):
                raise ValueError(
                    f"{location}: goto does not reach terminal cleanup before another optimizer call"
                )
            return_seen = True
            continue
        if return_seen:
            raise ValueError(f"{location}: executable statement follows return inside guard")
        if not any(pattern.match(line) for pattern in SAFE_CLEANUP_RES):
            raise ValueError(
                f"{location}: failure branch contains non-cleanup work: {line.strip()}"
            )

    if end_index < 0 or depth != 0:
        raise ValueError(f"{location}: unterminated failure guard")
    if not return_seen:
        raise ValueError(f"{location}: failure guard does not return")


def scan_file(path: Path, source_root: Path) -> list[tuple[str, str, int, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    current_subroutine = ""
    interface_depth = 0
    coalesced_calls: set[int] = set()
    found: list[tuple[str, str, int, str]] = []
    for index, line in enumerate(lines):
        if not significant(line):
            continue
        if END_INTERFACE_RE.match(line):
            interface_depth -= 1
            if interface_depth < 0:
                raise ValueError(f"{path}:{index + 1}: unmatched end interface")
            continue
        if INTERFACE_RE.match(line):
            interface_depth += 1
            continue
        match = SUBROUTINE_RE.match(line)
        if match and interface_depth == 0:
            current_subroutine = match.group(1).lower()
        if interface_depth:
            continue
        inline_refine = bool(INLINE_REFINE_RE.match(line))
        kind = (
            "direct" if DIRECT_RE.match(line) else
            "refine" if REFINE_RE.match(line) or inline_refine else ""
        )
        if not kind:
            continue
        if index in coalesced_calls:
            continue
        if not current_subroutine:
            raise ValueError(f"{path}:{index + 1}: call is outside a named subroutine")
        relative = path.relative_to(source_root).as_posix()
        location = f"{relative}:{index + 1} ({current_subroutine})"
        guard_index, alternate_calls = guard_after_call(
            lines, index, kind, inline_refine, location
        )
        validate_guard(lines, guard_index, location)
        coalesced_calls.update(alternate_calls)
        found.append((kind, current_subroutine, index + 1, relative))
    if interface_depth:
        raise ValueError(f"{path}: unterminated interface block")
    return found


def main() -> int:
    if len(sys.argv) != 2:
        print(f"usage: {sys.argv[0]} SOURCE_ROOT", file=sys.stderr)
        return 64
    source_root = Path(sys.argv[1]).resolve()
    source_dir = source_root / "src"
    if not source_dir.is_dir():
        print(f"missing source directory: {source_dir}", file=sys.stderr)
        return 65

    found: list[tuple[str, str, int, str]] = []
    try:
        for path in sorted(source_dir.rglob("*")):
            if path.is_file() and path.suffix.lower() in {".f", ".f90"}:
                found.extend(scan_file(path, source_root))
    except (OSError, UnicodeError, ValueError) as exc:
        print(f"OPTIMIZER_AFFINITY_CALLER_ORACLE_FAIL: {exc}", file=sys.stderr)
        return 1

    direct_counts = Counter(name for kind, name, _, _ in found if kind == "direct")
    refine_counts = Counter(name for kind, name, _, _ in found if kind == "refine")
    if direct_counts != DIRECT_EXPECTED:
        print(
            "OPTIMIZER_AFFINITY_CALLER_ORACLE_FAIL: direct caller inventory mismatch: "
            f"observed={dict(direct_counts)} expected={dict(DIRECT_EXPECTED)}",
            file=sys.stderr,
        )
        return 1
    if refine_counts != REFINE_EXPECTED:
        print(
            "OPTIMIZER_AFFINITY_CALLER_ORACLE_FAIL: refine propagation inventory mismatch: "
            f"observed={dict(refine_counts)} expected={dict(REFINE_EXPECTED)}",
            file=sys.stderr,
        )
        return 1

    ordinals: Counter[tuple[str, str]] = Counter()
    for kind, name, line, relative in found:
        ordinals[(kind, name)] += 1
        print(
            "CALLER_FAILURE_ORACLE_PASS "
            f"kind={kind} caller={name} occurrence={ordinals[(kind, name)]} "
            f"location={relative}:{line} post_failure_work=none"
        )
    print(
        "OPTIMIZER_AFFINITY_CALLER_ORACLE_PASS "
        f"direct={sum(direct_counts.values())} refine_propagation={sum(refine_counts.values())}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
