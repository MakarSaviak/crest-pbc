"""
fmlip_relay.backends._template_ase_backend
-------------------------------------------
Template for implementing a new ASE-calculator-based backend.

To add a real backend:
  1. Copy this file and rename it (e.g. ``my_calculator.py``).
  2. Replace ``MyCalculator`` / ``MyBackend`` with your names throughout.
  3. Implement ``__init__``: construct ``self._calc`` with your ASE calculator.
  4. Implement the ``name`` property.
  5. Register your backend in ``backends/__init__.py`` following the existing
     pattern (preferably as a lazy import so the dependency is optional).

Install requirements:
    pip install ".[my-extra]"   # adjust as needed

Usage (via CLI):
    fmlip-relay-server --port 54321 --backend my_backend \\
        --my-option value
"""

from __future__ import annotations

from ._ase_base import _ASEComputeMixin


class MyBackend(_ASEComputeMixin):
    """
    Wraps ``my_package.MyCalculator`` as an fmlip-relay backend.

    Parameters
    ----------
    model_path : str
        Path to the model file (adjust or remove as needed).
    device : str
        Device string, e.g. ``"cpu"`` or ``"cuda"``.
    """

    def __init__(self,
                 model_path: str,
                 device: str = "cpu"):
        try:
            from my_package import MyCalculator
        except ImportError as exc:
            raise ImportError(
                "my_package is not installed. Run: pip install my-package"
            ) from exc

        # ── construct the ASE calculator and store it as self._calc ──────────
        self._calc = MyCalculator(
            model=model_path,
            device=device,
        )
        self._model_path = model_path

    @property
    def name(self) -> str:
        return f"my_backend({self._model_path})"
