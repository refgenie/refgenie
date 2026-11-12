"""Refgenie: reference genome asset manager.

Top-level names are resolved lazily (PEP 562), so ``import refgenie`` does
not pay about 0.6s to load SQLModel, alembic, refget and cryptography.
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from .core import Refgenie
    from .managers.sources.manager import IndexFile
    from .models import BuildParams, GenomeAlias, GenomeDigest

__all__ = [
    "Refgenie",
    "BuildParams",
    "GenomeAlias",
    "GenomeDigest",
    "IndexFile",
]

# name -> (submodule, attribute)
_LAZY_ATTRS = {
    "Refgenie": ("refgenie.core", "Refgenie"),
    "BuildParams": ("refgenie.models", "BuildParams"),
    "GenomeAlias": ("refgenie.models", "GenomeAlias"),
    "GenomeDigest": ("refgenie.models", "GenomeDigest"),
    "IndexFile": ("refgenie.managers.sources.manager", "IndexFile"),
}


def _read_version() -> str:
    """Read the installed distribution's version, or "unknown" if not installed."""
    from importlib.metadata import PackageNotFoundError, version

    try:
        return version("refgenie")
    except PackageNotFoundError:
        return "unknown"


def __getattr__(name: str):
    if name == "__version__":
        value = _read_version()
        globals()[name] = value
        return value
    try:
        module_name, attr_name = _LAZY_ATTRS[name]
    except KeyError:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from None
    from importlib import import_module

    value = getattr(import_module(module_name), attr_name)
    globals()[name] = value  # cache so subsequent lookups skip __getattr__
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__) | {"__version__"})
