"""
Shared FastAPI dependencies for the refgenie dashboard and server.

Provides a singleton Refgenie instance and database session dependency
that all routers use via FastAPI's Depends() mechanism.
"""

from collections.abc import Generator
from typing import Annotated

from fastapi import Depends, Request
from fastapi.exceptions import RequestValidationError
from pydantic import Field
from sqlmodel import Session

from refgenie.const import SQL_INT_MAX
from refgenie.core import Refgenie

#: A database row id taken from a path or query parameter. Ids are PostgreSQL
#: INTEGER columns: a number outside that range used to reach the database and
#: come back as a 500. With this type FastAPI answers 422 first.
DbId = Annotated[int, Field(ge=0, le=SQL_INT_MAX)]

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


def reject_nul_bytes(request: Request) -> None:
    """422 for any path or query parameter that holds a NUL byte.

    Installed on the whole app (``create_app``), so it covers every route.
    PostgreSQL text cannot hold a NUL, so such a value used to reach the driver
    and come back as a 500. Raising ``RequestValidationError`` sends it through
    the same 422 handler as any other bad parameter.
    """
    params = [("path", k, v) for k, v in request.path_params.items()]
    params += [("query", k, v) for k, v in request.query_params.multi_items()]
    for where, key, value in params:
        if isinstance(value, str) and "\x00" in value:
            raise RequestValidationError(
                [
                    {
                        "type": "value_error",
                        "loc": (where, key),
                        "msg": "Value must not contain a NUL byte",
                        "input": value,
                    }
                ]
            )
