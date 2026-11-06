"""GA4GH page-based paging for the HTTP list routes.

Routes take ``page`` (from 0) and ``page_size`` and answer
``{"results": [...], "pagination": {"page", "page_size", "total"}}``. A
``page_size`` over the cap, or a page past the last one, is a 400.
"""

from collections.abc import Callable
from dataclasses import dataclass
from typing import Any, Generic, TypeVar

from fastapi import HTTPException, Query
from pydantic import BaseModel, Field

from refgenie.const import DEFAULT_PAGE_SIZE, MAX_PAGE_SIZE, SQL_BIGINT_MAX
from refgenie.models import Page

T = TypeVar("T")


class PageMeta(BaseModel):
    """Where a page sits in its list."""

    page: int = Field(description="This page's number, from 0")
    page_size: int = Field(description="Most results per page")
    total: int = Field(description="Results in the whole list")


class PageResponse(BaseModel, Generic[T]):
    """One page of a list."""

    results: list[T] = Field(description="The results on this page")
    pagination: PageMeta


@dataclass(frozen=True)
class PageParams:
    page: int
    page_size: int

    @property
    def offset(self) -> int:
        return self.page * self.page_size


def get_page_params(
    page: int = Query(
        0, ge=0, le=SQL_BIGINT_MAX // MAX_PAGE_SIZE, description="The page to return, from 0"
    ),
    page_size: int = Query(
        DEFAULT_PAGE_SIZE, ge=1, description=f"Results per page, at most {MAX_PAGE_SIZE}"
    ),
) -> PageParams:
    """FastAPI dependency for the paging query parameters.

    Raises:
        HTTPException: 400 if ``page_size`` is over the cap.
    """
    if page_size > MAX_PAGE_SIZE:
        raise HTTPException(400, f"page_size may be at most {MAX_PAGE_SIZE}")
    return PageParams(page=page, page_size=page_size)


def respond(page: Page, params: PageParams, convert: Callable[[Any], T]) -> PageResponse[T]:
    """The envelope for a manager's page, each item turned into its wire shape.

    Raises:
        HTTPException: 400 if ``params.page`` is past the last page. Page 0 of
            an empty list is not.
    """
    if params.page > 0 and params.offset >= page.total:
        last = max(0, (page.total - 1) // params.page_size)
        raise HTTPException(400, f"page {params.page} is past the last page ({last})")
    return PageResponse(
        results=[convert(item) for item in page.items],
        pagination=PageMeta(page=params.page, page_size=params.page_size, total=page.total),
    )
