"""
fmlip_relay.backends
-------------------
Backend registry.  Add new backends here so the CLI can find them by name.
"""

from .base  import BackendBase
from .dummy import DummyBackend

# Lazily registered so that missing optional deps don't break the import
_REGISTRY: dict[str, type[BackendBase]] = {
    "dummy": DummyBackend,
}

def _register_optional() -> None:
    from .lj import LJBackend
    _REGISTRY["lj"] = LJBackend

    try:
        from .mace import MACEBackend
        _REGISTRY["mace"] = MACEBackend
    except ImportError:
        pass

    try:
        from .mace_mp import MACEMPBackend
        _REGISTRY["mace_mp"] = MACEMPBackend
    except ImportError:
        pass

    try:
        from .mace_off import MACEOFFBackend
        _REGISTRY["mace_off"] = MACEOFFBackend
    except ImportError:
        pass

    try:
        from .uma import UMABackend
        _REGISTRY["uma"] = UMABackend
    except ImportError:
        pass

_register_optional()


def get_backend_class(name: str) -> type[BackendBase]:
    """Return the backend class for *name*, raising ``KeyError`` if unknown."""
    cls = _REGISTRY.get(name.lower())
    if cls is None:
        available = ", ".join(sorted(_REGISTRY))
        raise KeyError(
            f"Unknown backend '{name}'. Available backends: {available}"
        )
    return cls


__all__ = ["BackendBase", "DummyBackend", "get_backend_class"]
