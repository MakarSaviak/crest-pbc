"""
fmlip_relay.backends.mace_mp
---------------------------
Backend for the MACE-MP family of foundation models trained on the
Materials Project MPTrj dataset.  Covers 89 elements.

Models are downloaded automatically on first use and cached in
``~/.cache/mace`` (override with ``$XDG_CACHE_HOME``).

Install requirements:
    pip install ".[mace]"

Usage (via CLI):
    fmlip-relay-server --port 54321 --backend mace_mp
    fmlip-relay-server --port 54321 --backend mace_mp \\
        --mace-model medium-mpa-0 \\
        [--device cpu|cuda|cuda:0] \\
        [--dtype float32|float64] \\
        [--dispersion]

Available model sizes (as of mace-torch 0.3.10+):
    small, medium, large
    small-0b, medium-0b, small-0b2, medium-0b2, medium-0b3
    medium-mpa-0  (default)
    large-0b2
    medium-omat-0

See https://github.com/ACEsuit/mace for the latest releases.
"""

from __future__ import annotations

from ._ase_base import _ASEComputeMixin

# Kept here so the CLI can show them in help text without importing torch
_KNOWN_MODELS = (
    "small", "medium", "large",
    "small-0b", "medium-0b",
    "small-0b2", "medium-0b2", "medium-0b3",
    "medium-mpa-0",
    "large-0b2",
    "medium-omat-0",
)
_DEFAULT_MODEL = "medium-mpa-0"


class MACEMPBackend(_ASEComputeMixin):
    """
    Wraps ``mace.calculators.mace_mp`` foundation model calculator.

    Parameters
    ----------
    model : str
        Model size or release tag, e.g. ``"medium-mpa-0"`` (default),
        ``"small"``, ``"large"``.  See module docstring for all options.
    device : str
        PyTorch device string.
    dtype : str
        ``"float32"`` (default for foundation models) or ``"float64"``.
    dispersion : bool
        Add a D3 dispersion correction on top of the MACE energy.
        Requires ``torch-dftd`` to be installed.
    """

    def __init__(self,
                 model:      str  = _DEFAULT_MODEL,
                 device:     str  = "cpu",
                 dtype:      str  = "float32",
                 dispersion: bool = False):
        try:
            from mace.calculators import mace_mp
        except ImportError as exc:
            raise ImportError(
                "MACE is not installed. Run: pip install mace-torch"
            ) from exc

        self._calc       = mace_mp(
            model=model,
            device=device,
            default_dtype=dtype,
            dispersion=dispersion,
        )
        self._model      = model
        self._dispersion = dispersion
        self._device     = device
        self._dtype      = dtype

    @property
    def name(self) -> str:
        disp = "+D3" if self._dispersion else ""
        return f"mace_mp({self._model}{disp})"
