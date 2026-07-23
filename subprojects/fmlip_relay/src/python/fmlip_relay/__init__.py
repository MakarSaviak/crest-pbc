"""
fmlip_relay
----------
Persistent ML potential server for Fortran clients.
"""

__version__ = "0.1.0"

from .server   import run
from .backends import get_backend_class, BackendBase

__all__ = ["run", "get_backend_class", "BackendBase", "__version__"]
