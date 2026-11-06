"""Shared FastAPI dependencies: the Refgenie instance and parameter guards."""

from typing import Annotated

from fastapi import Request
from fastapi.exceptions import RequestValidationError
from pydantic import Field

from refgenie.const import SQL_INT_MAX
from refgenie.core import Refgenie

#: A database row id from a path or query parameter. Out-of-range values are a
#: 422, not a database error.
DbId = Annotated[int, Field(ge=0, le=SQL_INT_MAX)]

_refgenie_instance: Refgenie | None = None


def get_refgenie() -> Refgenie:
    """The Refgenie instance shared by all routes."""
    global _refgenie_instance
    if _refgenie_instance is None:
        _refgenie_instance = Refgenie()
    return _refgenie_instance


def reject_nul_bytes(request: Request) -> None:
    """422 for any path or query parameter that holds a NUL byte.

    Installed app-wide by ``create_app``. PostgreSQL text cannot hold a NUL.
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
