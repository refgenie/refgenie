"""Refgenie: reference genome asset manager.

The public top-level names are resolved lazily (PEP 562). Importing them
eagerly pulled in SQLModel, alembic, refget and cryptography on every
``import refgenie`` -- about 0.6s -- which every CLI invocation and every
test subprocess paid even when it only needed ``refgenie.const``.

``from refgenie import Refgenie`` still works; it just does the expensive
work at attribute-access time instead of at package-import time.
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover - import-time cost is the whole point
    from .managers.sources.manager import IndexFile
    from .models import BuildParams
    from .populator import looper_refgenie_populate_local
    from .core import Refgenie

__all__ = [
    "Refgenie",
    "BuildParams",
    "IndexFile",
    "looper_refgenie_populate_local",
]

# name -> (submodule, attribute)
_LAZY_ATTRS = {
    "Refgenie": ("refgenie.core", "Refgenie"),
    "BuildParams": ("refgenie.models", "BuildParams"),
    "IndexFile": ("refgenie.managers.sources.manager", "IndexFile"),
    "looper_refgenie_populate_local": (
        "refgenie.populator",
        "looper_refgenie_populate_local",
    ),
}


def _read_version() -> str:
    """Read the installed distribution's version.

    Kept lazy for the same reason as everything else here: ``importlib.metadata``
    walks ``sys.path`` to find the distribution, and nothing on the CLI's hot
    path needs the version. Reading it from metadata rather than hardcoding a
    literal keeps pyproject.toml the single source of truth.
    """
    from importlib.metadata import PackageNotFoundError, version

    try:
        return version("refgenie")
    except PackageNotFoundError:
        # Imported from a source tree that was never installed (e.g. a bare
        # `python -c` from the repo root). Not an error worth raising for.
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
