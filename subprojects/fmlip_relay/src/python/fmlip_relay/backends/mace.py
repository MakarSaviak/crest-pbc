"""
fmlip_relay.backends.mace
------------------------
Backend for custom-trained MACE models loaded from a local file.

Install requirements:
    pip install ".[mace]"

Usage (via CLI):
    fmlip-relay-server --port 54321 --backend mace \\
        --model /path/to/model.model \\
        [--device cpu|cuda|cuda:0] \\
        [--dtype float32|float64]
"""

from __future__ import annotations

from ._ase_base import _ASEComputeMixin


class MACEBackend(_ASEComputeMixin):
    """
    Wraps ``mace.calculators.MACECalculator`` for a custom-trained model.

    Parameters
    ----------
    model_path : str
        Path to the ``.model`` file produced by MACE training.
    device : str
        PyTorch device string, e.g. ``"cpu"``, ``"cuda"``, ``"cuda:0"``.
    dtype : str
        ``"float64"`` (default) or ``"float32"``.
    """

    def __init__(self,
                 model_path: str,
                 device: str = "cpu",
                 dtype:  str = "float64"):
        try:
            from mace.calculators import MACECalculator
        except ImportError as exc:
            raise ImportError(
                "MACE is not installed. Run: pip install mace-torch"
            ) from exc

        self._calc       = MACECalculator(
            model_paths=model_path,
            device=device,
            default_dtype=dtype,
        )
        self._model_path = model_path
        self._device     = device
        self._dtype      = dtype

    @property
    def name(self) -> str:
        return f"mace({self._model_path})"
