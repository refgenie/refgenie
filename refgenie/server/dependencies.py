"""
Shared FastAPI dependencies for the refgenie dashboard and server.

Provides a singleton Refgenie instance and database session dependency
that all routers use via FastAPI's Depends() mechanism.
"""

from collections.abc import Generator

from fastapi import Depends
from sqlmodel import Session

from refgenie.core import Refgenie

_refgenie_instance: Refgenie | None = None


def get_refgenie() -> Refgenie:
    """Singleton Refgenie instance shared across all dashboard endpoints."""
    global _refgenie_instance
    if _refgenie_instance is None:
        _refgenie_instance = Refgenie()
    return _refgenie_instance


def get_db_session(rgc: Refgenie = Depends(get_refgenie)) -> Generator[Session, None, None]:
    """FastAPI dependency that yields a database session from the shared Refgenie instance.

    ``get_refgenie`` must arrive through ``Depends`` rather than be called
    directly: a plain call is invisible to ``app.dependency_overrides``, so a
    handler taking both ``rgc`` and ``session`` would be handed the overridden
    Refgenie alongside a session bound to the default one -- two databases.
    """
    with Session(rgc.database_engine) as session:
        yield session
