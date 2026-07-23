"""
fmlip_relay.backends.uma
------------------------
Backend for Meta FAIR's UMA (Universal Models for Atoms) foundation models,
provided through the ``fairchem`` (v2) package.

UMA is a single multi-task model covering molecules, materials, catalysis,
molecular crystals, and MOFs.  The active domain is selected at construction
time via ``task_name``:

    omol   organic / inorganic molecules (uses charge & spin)
    omat   inorganic bulk materials
    omc    molecular crystals
    oc20   heterogeneous catalysis (adsorbates on surfaces)
    odac   metal-organic frameworks / direct-air capture

For the ``omol`` task the total charge and spin multiplicity are read from
``atoms.info["charge"]`` / ``atoms.info["spin"]`` — these are forwarded from
the relay protocol automatically by the shared ASE mixin.

The UMA checkpoints are gated on the Hugging Face Hub.  Before first use you
must request access on the model page and authenticate locally, e.g.:

    huggingface-cli login          # or set $HF_TOKEN

Models are then downloaded automatically and cached in
``~/.cache/huggingface`` (override with ``$HF_HOME``).

Install requirements:
    pip install ".[uma]"

Usage (via CLI):
    fmlip-relay-server --port 54321 --backend uma --uma-task omol
    fmlip-relay-server --port 54321 --backend uma \\
        --uma-model uma-s-1p2 \\
        --uma-task omat \\
        [--device cpu|cuda|cuda:0]

Available model variants (as of fairchem-core 2.x):
    uma-s-1, uma-s-1p1, uma-s-1p2 (default), uma-m-1

See https://github.com/facebookresearch/fairchem for the latest releases.
"""

from __future__ import annotations

from ._ase_base import _ASEComputeMixin

# Kept here so the CLI can show them in help text without importing torch
_KNOWN_MODELS = ("uma-s-1", "uma-s-1p1", "uma-s-1p2", "uma-m-1")
_DEFAULT_MODEL = "uma-s-1p2"

_KNOWN_TASKS = ("omol", "omat", "omc", "oc20", "odac")
_DEFAULT_TASK = "omol"


class UMABackend(_ASEComputeMixin):
    """
    Wraps ``fairchem.core.FAIRChemCalculator`` for a UMA foundation model.

    Parameters
    ----------
    model : str
        UMA checkpoint name, e.g. ``"uma-s-1p2"`` (default), ``"uma-s-1p1"``,
        ``"uma-s-1"``, ``"uma-m-1"``.  Resolved against the Hugging Face Hub.
    task : str
        Task / domain head to evaluate: one of ``omol``, ``omat``, ``omc``,
        ``oc20``, ``odac``.  Defaults to ``"omol"`` (molecules; honours the
        charge and spin passed through the protocol).
    device : str
        PyTorch device string, e.g. ``"cpu"``, ``"cuda"``, ``"cuda:0"``.
        UMA runs best on a GPU; CPU evaluation works but is slow.
    """

    def __init__(self,
                 model:  str = _DEFAULT_MODEL,
                 task:   str = _DEFAULT_TASK,
                 device: str = "cpu"):
        try:
            from fairchem.core import pretrained_mlip, FAIRChemCalculator
        except ImportError as exc:
            raise ImportError(
                "fairchem is not installed. Run: pip install fairchem-core"
            ) from exc

        if task not in _KNOWN_TASKS:
            raise ValueError(
                f"Unknown uma task '{task}'. "
                f"Choose from: {', '.join(_KNOWN_TASKS)}"
            )

        # ── load the shared multi-task predict unit, then bind a task head ───
        predictor = pretrained_mlip.get_predict_unit(model, device=device)
        self._calc = FAIRChemCalculator(predictor, task_name=task)

        self._model  = model
        self._task   = task
        self._device = device
        self._dtype  = "float32"   # fairchem runs UMA inference in float32

    @property
    def name(self) -> str:
        return f"uma({self._model}/{self._task})"
