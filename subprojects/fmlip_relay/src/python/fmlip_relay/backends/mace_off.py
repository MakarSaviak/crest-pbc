"""
fmlip_relay.backends.mace_off
----------------------------
Backend for the MACE-OFF23 family of organic force fields.

Parameterised for 10 elements: H, C, N, O, P, S, F, Cl, Br, I.
Suitable for neutral organic molecules in gas phase, liquid phase, or
organic crystals.  Published under the Academic Software License (ASL).

Models are downloaded automatically on first use and cached in
``~/.cache/mace`` (override with ``$XDG_CACHE_HOME``).

Install requirements:
    pip install ".[mace]"

Usage (via CLI):
    fmlip-relay-server --port 54321 --backend mace_off
    fmlip-relay-server --port 54321 --backend mace_off \\
        --mace-model medium \\
        [--device cpu|cuda|cuda:0] \\
        [--dtype float32|float64]

Available model sizes:
    small, medium (default), large
"""

from __future__ import annotations

from ._ase_base import _ASEComputeMixin

_KNOWN_MODELS = ("small", "medium", "large")
_DEFAULT_MODEL = "medium"


class MACEOFFBackend(_ASEComputeMixin):
    """
    Wraps ``mace.calculators.mace_off`` organic foundation model calculator.

    Parameters
    ----------
    model : str
        Model size: ``"small"``, ``"medium"`` (default), or ``"large"``.
    device : str
        PyTorch device string.
    dtype : str
        ``"float32"`` (default) or ``"float64"``.
    """

    def __init__(self,
                 model:  str = _DEFAULT_MODEL,
                 device: str = "cpu",
                 dtype:  str = "float32"):
        try:
            from mace.calculators import mace_off
        except ImportError as exc:
            raise ImportError(
                "MACE is not installed. Run: pip install mace-torch"
            ) from exc

        if model not in _KNOWN_MODELS:
            raise ValueError(
                f"Unknown mace_off model '{model}'. "
                f"Choose from: {', '.join(_KNOWN_MODELS)}"
            )

        self._calc   = mace_off(
            model=model,
            device=device,
            default_dtype=dtype,
        )
        self._model  = model
        self._device = device
        self._dtype  = dtype

    @property
    def name(self) -> str:
        return f"mace_off({self._model})"
