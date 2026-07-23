"""
fmlip_relay.check
-----------------
Environment checker: probes every known backend and reports which ones
are usable in the current Python environment.

Invoked via:
    fmlip-relay-check
    python -m fmlip_relay check
"""

from __future__ import annotations

import importlib
import importlib.metadata
import sys


# ── colour helpers (no external deps) ────────────────────────────────────────

def _supports_colour() -> bool:
    return hasattr(sys.stdout, "isatty") and sys.stdout.isatty()


_COLOUR = _supports_colour()
_GREEN  = "\033[32m" if _COLOUR else ""
_RED    = "\033[31m" if _COLOUR else ""
_YELLOW = "\033[33m" if _COLOUR else ""
_BOLD   = "\033[1m"  if _COLOUR else ""
_RESET  = "\033[0m"  if _COLOUR else ""

OK   = f"{_GREEN}✓ ok{_RESET}"
FAIL = f"{_RED}✗ unavailable{_RESET}"
WARN = f"{_YELLOW}~ partial{_RESET}"


# ── individual probe functions ────────────────────────────────────────────────

def _probe_numpy() -> tuple[bool, str]:
    try:
        import numpy as np
        return True, f"numpy {np.__version__}"
    except ImportError as e:
        return False, str(e)


def _probe_torch() -> tuple[bool, str]:
    try:
        import torch
        cuda = f"  CUDA available: {torch.cuda.is_available()}"
        return True, f"torch {torch.__version__}\n    {cuda}"
    except ImportError as e:
        return False, str(e)


def _probe_mace() -> tuple[bool, str]:
    try:
        importlib.import_module("mace")
        ver = importlib.metadata.version("mace-torch")
        return True, f"mace-torch {ver}"
    except ImportError as e:
        return False, str(e)
    except importlib.metadata.PackageNotFoundError:
        return True, "mace (version unknown)"


def _probe_fairchem() -> tuple[bool, str]:
    try:
        importlib.import_module("fairchem.core")
        ver = importlib.metadata.version("fairchem-core")
        return True, f"fairchem-core {ver}"
    except ImportError as e:
        return False, str(e)
    except importlib.metadata.PackageNotFoundError:
        return True, "fairchem (version unknown)"


def _probe_ase() -> tuple[bool, str]:
    try:
        import ase
        return True, f"ase {ase.__version__}"
    except ImportError as e:
        return False, str(e)


def _probe_dispersion() -> tuple[bool, str]:
    try:
        importlib.import_module("torch_dftd")
        ver = importlib.metadata.version("torch-dftd")
        return True, f"torch-dftd {ver}"
    except ImportError as e:
        return False, str(e)
    except importlib.metadata.PackageNotFoundError:
        return True, "torch-dftd (version unknown)"


def _probe_backend_import(module: str, cls: str) -> tuple[bool, str]:
    """Try to import a backend class without instantiating it."""
    try:
        mod = importlib.import_module(module)
        getattr(mod, cls)
        return True, ""
    except Exception as e:
        return False, str(e)


# ── backend table ─────────────────────────────────────────────────────────────

# Each entry: (cli_name, module, class, description, extra_probe_fn | None)
_BACKENDS = [
    ("dummy",    "fmlip_relay.backends.dummy",   "DummyBackend",  "Random numbers, no deps",          None),
    ("lj",       "fmlip_relay.backends.lj",      "LJBackend",     "Lennard-Jones, no deps",            None),
    ("mace",     "fmlip_relay.backends.mace",    "MACEBackend",   "Custom MACE model from file",       None),
    ("mace_mp",  "fmlip_relay.backends.mace_mp", "MACEMPBackend", "MACE-MP foundation (89 elements)",  None),
    ("mace_off", "fmlip_relay.backends.mace_off","MACEOFFBackend","MACE-OFF foundation (organic)",      None),
    ("uma",      "fmlip_relay.backends.uma",     "UMABackend",    "UMA universal model (fairchem)",    None),
]


# ── report ────────────────────────────────────────────────────────────────────

def run_check() -> int:
    """
    Run all probes and print a formatted report.
    Returns 0 if all required deps are present, 1 otherwise.
    """
    print(f"\n{_BOLD}fmlip-relay environment check{_RESET}")
    print("=" * 52)

    # ── core dependencies ─────────────────────────────────────────────────────
    print(f"\n{_BOLD}Core dependencies{_RESET}")
    numpy_ok, numpy_info = _probe_numpy()
    _print_row("numpy", numpy_ok, numpy_info)

    # ── ML dependencies ───────────────────────────────────────────────────────
    print(f"\n{_BOLD}ML dependencies{_RESET}")
    torch_ok,      torch_info      = _probe_torch()
    mace_ok,       mace_info       = _probe_mace()
    fairchem_ok,   fairchem_info   = _probe_fairchem()
    ase_ok,        ase_info        = _probe_ase()
    dispersion_ok, dispersion_info = _probe_dispersion()

    _print_row("torch",      torch_ok,      torch_info)
    _print_row("mace-torch", mace_ok,       mace_info)
    _print_row("fairchem",   fairchem_ok,   fairchem_info,
               note="needed for the uma backend")
    _print_row("ase",        ase_ok,        ase_info)
    _print_row("torch-dftd", dispersion_ok, dispersion_info,
               note="optional — needed for --dispersion with mace_mp")

    # ── backends ──────────────────────────────────────────────────────────────
    print(f"\n{_BOLD}Available backends{_RESET}")
    all_ok = True
    for cli_name, module, cls, description, _ in _BACKENDS:
        ok, err = _probe_backend_import(module, cls)
        if not ok:
            all_ok = False
        _print_backend_row(cli_name, ok, description, err)

    # ── summary ───────────────────────────────────────────────────────────────
    print()
    if not numpy_ok:
        print(f"{_RED}ERROR: numpy is missing — nothing will work.{_RESET}")
        return 1
    if not (torch_ok and mace_ok and ase_ok):
        print(f"{_YELLOW}MACE backends unavailable. "
              f"Install with:  pip install \"fmlip-relay[mace]\"{_RESET}")
    if not (torch_ok and fairchem_ok and ase_ok):
        print(f"{_YELLOW}UMA backend unavailable. "
              f"Install with:  pip install \"fmlip-relay[uma]\"  "
              f"(then 'huggingface-cli login' for gated checkpoints){_RESET}")
    if not dispersion_ok:
        print(f"{_YELLOW}D3 dispersion unavailable. "
              f"Install with:  pip install torch-dftd{_RESET}")
    if all_ok or (numpy_ok):
        print(f"{_GREEN}Core is functional.{_RESET}")

    return 0


def _print_row(label: str, ok: bool, info: str, note: str = "") -> None:
    status = OK if ok else FAIL
    # indent continuation lines of info
    info_fmt = info.replace("\n", "\n      ") if info else ""
    print(f"  {status}  {label:<14} {info_fmt}")
    if note:
        print(f"             {_YELLOW}({note}){_RESET}")


def _print_backend_row(name: str, ok: bool, description: str, err: str) -> None:
    status = OK if ok else FAIL
    print(f"  {status}  --backend {name:<10}  {description}")
    if err:
        # truncate long import error chains to the last line
        last_line = err.strip().splitlines()[-1]
        print(f"             {_RED}{last_line}{_RESET}")


def main() -> None:
    sys.exit(run_check())


if __name__ == "__main__":
    main()
