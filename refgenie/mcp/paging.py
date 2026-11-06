"""
GA4GH token-based paging for the MCP list tools.

A tool takes ``page_size`` (default 20, at most 100; out-of-range values are
clamped, and the size used is echoed back) and an opaque ``token``. The token
is the urlsafe base64 of ``{"v": 1, "o": offset, "h": filters_hash}``: the
offset of the next page, and a hash of the tool and its filters, so a token
from one query is refused by another instead of paging through the wrong list.
"""

import base64
import binascii
import hashlib
import json
from typing import Generic, TypeVar

from mcp.server.mcpserver.exceptions import ToolError
from pydantic import BaseModel

from refgenie.models import Page

PAGE_SIZE_DEFAULT, PAGE_SIZE_CAP = 20, 100

R = TypeVar("R")


class TokenPagination(BaseModel):
    page_size: int
    next_page_token: str | None
    total: int


class ResultPage(BaseModel, Generic[R]):
    results: list[R]
    pagination: TokenPagination
    hint: str | None = None


def clamp(page_size: int | None) -> int:
    """``page_size`` forced into ``[1, PAGE_SIZE_CAP]``; None or 0 means the default."""
    return max(1, min(page_size or PAGE_SIZE_DEFAULT, PAGE_SIZE_CAP))


def filters_hash(tool: str, **filters) -> str:
    """10 hex characters naming one tool and its non-None filters."""
    kept = {k: v for k, v in sorted(filters.items()) if v is not None}
    blob = json.dumps([tool, kept], sort_keys=True, default=str)
    return hashlib.sha256(blob.encode()).hexdigest()[:10]


def encode(offset: int, fh: str) -> str:
    """The token for the page starting at ``offset`` of the query hashed ``fh``."""
    raw = json.dumps({"v": 1, "o": offset, "h": fh}, separators=(",", ":")).encode()
    return base64.urlsafe_b64encode(raw).decode().rstrip("=")


def decode(token: str | None, fh: str) -> int:
    """
    The offset a token points at; 0 for no token.

    Raises:
        ToolError: If the token was not issued here, or was issued for other filters.
    """
    if not token:
        return 0
    try:
        payload = json.loads(base64.urlsafe_b64decode(token + "=" * (-len(token) % 4)))
        offset, issued_for = int(payload["o"]), payload["h"]
        if payload.get("v") != 1 or offset < 0:
            raise ValueError
    except (binascii.Error, ValueError, TypeError, KeyError, AttributeError):
        raise ToolError("token is not one this server issued") from None
    if issued_for != fh:
        raise ToolError("token does not match this query; start again without a token")
    return offset


def result_page(page: Page, rows: list, fh: str, hint: str | None = None) -> ResultPage:
    """The tool's envelope for one manager ``Page`` and its rows already converted."""
    end = page.offset + len(page.items)
    return ResultPage(
        results=rows,
        pagination=TokenPagination(
            page_size=page.limit,
            next_page_token=encode(end, fh) if end < page.total else None,
            total=page.total,
        ),
        hint=hint,
    )
