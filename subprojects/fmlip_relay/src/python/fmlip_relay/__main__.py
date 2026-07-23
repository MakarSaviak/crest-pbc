"""
fmlip_relay.__main__
-------------------
CLI entry point.  Installed as the ``fmlip-relay-server`` console script.
Also invokable as ``python -m fmlip_relay``.

Usage
-----
    # Custom-trained MACE model
    fmlip-relay-server --port 54321 --backend mace --model /path/to/model.model

    # MACE-MP foundation model (Materials Project, 89 elements)
    fmlip-relay-server --port 54321 --backend mace_mp
    fmlip-relay-server --port 54321 --backend mace_mp --mace-model large-0b2

    # MACE-OFF foundation model (organic, 10 elements)
    fmlip-relay-server --port 54321 --backend mace_off
    fmlip-relay-server --port 54321 --backend mace_off --mace-model small

    # UMA universal model (Meta FAIR / fairchem) — pick a task head
    fmlip-relay-server --port 54321 --backend uma --uma-task omol
    fmlip-relay-server --port 54321 --backend uma --uma-model uma-m-1 --uma-task omat

    # Lennard-Jones
    fmlip-relay-server --port 54321 --backend lj --lj-epsilon 0.0104

    # Dummy (no model required)
    fmlip-relay-server --port 54321 --backend dummy

    # Single-point test — loads geometry, runs one evaluation, prints summary, exits
    fmlip-relay-server --backend mace_mp --test molecule.xyz
    fmlip-relay-server --backend lj      --test crystal.xyz
"""

from __future__ import annotations

import argparse
import logging
import sys

from .backends import get_backend_class
from .backends.mace_mp  import _KNOWN_MODELS as _MACE_MP_MODELS,  _DEFAULT_MODEL as _MACE_MP_DEFAULT
from .backends.mace_off import _KNOWN_MODELS as _MACE_OFF_MODELS, _DEFAULT_MODEL as _MACE_OFF_DEFAULT
from .backends.uma      import (_KNOWN_MODELS as _UMA_MODELS, _DEFAULT_MODEL as _UMA_DEFAULT,
                                _KNOWN_TASKS  as _UMA_TASKS,  _DEFAULT_TASK  as _UMA_TASK_DEFAULT)
from . import server
from . import test_mode


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fmlip-relay-server",
        description="Persistent ML potential server for Fortran clients.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # ── required ──────────────────────────────────────────────────────────────
    p.add_argument("--port",    type=int, default=None,
                   help="TCP port to listen on (loopback only). "
                        "Required unless --test is used.")
    p.add_argument("--backend", type=str, required=True,
                   help="Backend name: mace | mace_mp | mace_off | uma | lj | dummy")

    # ── test mode ─────────────────────────────────────────────────────────────
    grp_test = p.add_argument_group("Single-point test options  (--test)")
    grp_test.add_argument("--test", type=str, default=None, metavar="GEOMETRY",
                          help="Run a single-point calculation on the given geometry "
                               "file (XYZ / extended-XYZ), print a summary, and exit. "
                               "--port is not required in this mode.")
    grp_test.add_argument("--charge", type=int, default=None, metavar="Q",
                          help="Total charge for the --test evaluation. "
                               "Default: value from the geometry file if present, "
                               "else 0 (neutral). Used by charge-aware backends "
                               "(e.g. UMA omol).")
    grp_test.add_argument("--spin", type=int, default=None, metavar="M",
                          help="Spin multiplicity 2S+1 for the --test evaluation. "
                               "Default: value from the geometry file if present, "
                               "else 1 (ground-state singlet). Used by spin-aware "
                               "backends (e.g. UMA omol).")

    # ── logging ───────────────────────────────────────────────────────────────
    p.add_argument("--log",      type=str, default=None,
                   help="Optional log file path (stderr is always used)")
    p.add_argument("--loglevel", type=str, default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"],
                   help="Log verbosity")

    # ── backend: mace (custom model file) ────────────────────────────────────
    grp_mace = p.add_argument_group(
        "Custom MACE model options  (--backend mace)"
    )
    grp_mace.add_argument("--model", type=str, default=None,
                          help="Path to .model file (required for mace backend)")

    # ── shared neural-network options ─────────────────────────────────────────
    grp_nn = p.add_argument_group(
        "Shared neural-network options  (mace | mace_mp | mace_off | uma)"
    )
    grp_nn.add_argument("--device", type=str, default="cpu",
                        help="PyTorch device string, e.g. cpu | cuda | cuda:0")

    # ── shared MACE precision option ──────────────────────────────────────────
    grp_mace_shared = p.add_argument_group(
        "Shared MACE options  (mace | mace_mp | mace_off)"
    )
    grp_mace_shared.add_argument("--dtype",  type=str, default=None,
                                 choices=["float32", "float64"],
                                 help="Floating-point precision "
                                      "(default: float64 for mace, float32 for foundation models). "
                                      "UMA precision is fixed by fairchem (float32).")

    # ── backend: mace_mp ─────────────────────────────────────────────────────
    grp_mp = p.add_argument_group(
        "MACE-MP options  (--backend mace_mp)  — Materials Project, 89 elements"
    )
    grp_mp.add_argument(
        "--mace-model", type=str, default=None,
        metavar="MODEL",
        help=(f"Foundation model variant. "
              f"mace_mp default: '{_MACE_MP_DEFAULT}', options: "
              f"{', '.join(_MACE_MP_MODELS)}. "
              f"mace_off default: '{_MACE_OFF_DEFAULT}', options: "
              f"{', '.join(_MACE_OFF_MODELS)}.")
    )
    grp_mp.add_argument("--dispersion", action="store_true", default=False,
                        help="Add D3 dispersion correction (mace_mp only; "
                             "requires torch-dftd)")

    # ── backend: uma ─────────────────────────────────────────────────────────
    grp_uma = p.add_argument_group(
        "UMA options  (--backend uma)  — Meta FAIR universal model (fairchem)"
    )
    grp_uma.add_argument(
        "--uma-model", type=str, default=None, metavar="MODEL",
        help=(f"UMA checkpoint. Default: '{_UMA_DEFAULT}', options: "
              f"{', '.join(_UMA_MODELS)}.")
    )
    grp_uma.add_argument(
        "--uma-task", type=str, default=None, metavar="TASK",
        choices=list(_UMA_TASKS),
        help=(f"Task / domain head. Default: '{_UMA_TASK_DEFAULT}', options: "
              f"{', '.join(_UMA_TASKS)}.")
    )

    # ── backend: lj ──────────────────────────────────────────────────────────
    grp_lj = p.add_argument_group("Lennard-Jones options  (--backend lj)")
    grp_lj.add_argument("--lj-epsilon", type=float, default=0.0104,
                        help="Well depth ε in eV")
    grp_lj.add_argument("--lj-sigma",   type=float, default=3.40,
                        help="Zero-crossing σ in Å")
    grp_lj.add_argument("--lj-cutoff",  type=float, default=8.50,
                        help="Cutoff radius in Å")

    return p


def _configure_logging(loglevel: str, logfile: str | None) -> None:
    handlers: list[logging.Handler] = [logging.StreamHandler(sys.stderr)]
    if logfile:
        handlers.append(logging.FileHandler(logfile))
    # Keep the root logger at WARNING so chatty third-party libraries (httpx,
    # fairchem, urllib3, huggingface_hub, …) don't spam at INFO during model
    # download/load.  Only our own package logs at the requested level; its
    # records still reach the root handlers via propagation.
    logging.basicConfig(
        level=logging.WARNING,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        handlers=handlers,
    )
    logging.getLogger("fmlip_relay").setLevel(getattr(logging, loglevel))


def _build_backend(args: argparse.Namespace):
    cls     = get_backend_class(args.backend)
    backend = args.backend.lower()
    kwargs: dict = {}

    if backend == "mace":
        if not args.model:
            raise ValueError("--model is required for the mace backend")
        kwargs = dict(
            model_path = args.model,
            device     = args.device,
            dtype      = args.dtype or "float64",
        )

    elif backend == "mace_mp":
        kwargs = dict(
            model      = args.mace_model or _MACE_MP_DEFAULT,
            device     = args.device,
            dtype      = args.dtype or "float32",
            dispersion = args.dispersion,
        )

    elif backend == "mace_off":
        if args.dispersion:
            raise ValueError("--dispersion is not supported for mace_off")
        kwargs = dict(
            model  = args.mace_model or _MACE_OFF_DEFAULT,
            device = args.device,
            dtype  = args.dtype or "float32",
        )

    elif backend == "uma":
        kwargs = dict(
            model  = args.uma_model or _UMA_DEFAULT,
            task   = args.uma_task  or _UMA_TASK_DEFAULT,
            device = args.device,
        )

    elif backend == "lj":
        kwargs = dict(
            epsilon = args.lj_epsilon,
            sigma   = args.lj_sigma,
            r_cut   = args.lj_cutoff,
        )

    return cls(**kwargs)


def main(argv: list[str] | None = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)

    _configure_logging(args.loglevel, args.log)
    log = logging.getLogger(__name__)

    # ── validate port requirement ─────────────────────────────────────────────
    if args.test is None and args.port is None:
        parser.error("--port is required when not using --test")

    # ── build backend (shared by both modes) ──────────────────────────────────
    try:
        backend = _build_backend(args)
    except (KeyError, ValueError) as exc:
        parser.error(str(exc))

    log.info("Backend '%s' loaded successfully.", backend.name)

    # ── test mode ─────────────────────────────────────────────────────────────
    if args.test is not None:
        sys.exit(test_mode.run_test(backend, args.test,
                                    charge=args.charge, spin=args.spin))

    # ── server mode ───────────────────────────────────────────────────────────
    server.run(args.port, backend)


if __name__ == "__main__":
    main()
