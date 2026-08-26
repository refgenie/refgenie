from typing import Any, Generic, TypeVar

from fastapi import Query
from pydantic import BaseModel, Field
from sqlalchemy import Select
from sqlmodel import Session, SQLModel, func, select

from refgenie.const import DEFAULT_PAGE_SIZE, MAX_PAGE_SIZE

T = TypeVar("T", bound=SQLModel)

_T = TypeVar("_T")


class PaginationMeta(BaseModel):
    """
    Metadata for paginated responses.
    """

    offset: int = Field(description="Starting index of the results")
    limit: int = Field(description="Maximum number of results per page")
    total: int = Field(description="Total number of items available")


class PaginatedResponse(BaseModel, Generic[_T]):
    """
    Generic paginated response model.
    """

    items: list[_T] = Field(description="List of items for this page")
    pagination: PaginationMeta = Field(description="Pagination metadata")


class PaginationParams(BaseModel):
    """
    Query parameters for pagination.
    """

    offset: int = Field(default=0, ge=0, description="Number of items to skip")
    limit: int = Field(
        default=DEFAULT_PAGE_SIZE,
        ge=1,
        le=MAX_PAGE_SIZE,
        description="Maximum number of items to return",
    )

    @property
    def skip(self) -> int:
        """Alias for offset to match SQLModel conventions."""
        return self.offset


def paginate_query(
    session: Session,
    query: Select[Any],
    pagination: PaginationParams,
) -> PaginatedResponse[Any]:
    """
    Apply pagination to a SQLModel query and return a paginated response.

    Args:
        session: Database session
        query: SQLModel select query
        pagination: Pagination parameters (offset, limit)

    Returns:
        PaginatedResponse with items and pagination metadata
    """
    count_query = select(func.count()).select_from(query.subquery())
    total = session.exec(count_query).one()

    paginated_query = query.offset(pagination.offset).limit(pagination.limit)
    items = session.exec(paginated_query).unique().all()  # type: ignore

    pagination_meta = PaginationMeta(
        offset=pagination.offset,
        limit=pagination.limit,
        total=total,
    )

    return PaginatedResponse(items=items, pagination=pagination_meta)


def paginate_list(
    items: list[T],
    pagination: PaginationParams,
) -> PaginatedResponse[T]:
    """
    Apply pagination to an in-memory list.

    Args:
        items: List of items to paginate
        pagination: Pagination parameters (offset, limit)

    Returns:
        PaginatedResponse with paginated items and metadata
    """
    total = len(items)
    start = pagination.offset
    end = start + pagination.limit

    paginated_items = items[start:end]

    pagination_meta = PaginationMeta(
        offset=pagination.offset,
        limit=pagination.limit,
        total=total,
    )

    return PaginatedResponse(items=paginated_items, pagination=pagination_meta)


def get_pagination_params(
    offset: int = Query(0, ge=0, description="Number of items to skip"),
    limit: int = Query(
        DEFAULT_PAGE_SIZE,
        ge=1,
        le=MAX_PAGE_SIZE,
        description="Maximum number of items to return",
    ),
) -> PaginationParams:
    """
    Dependency injection function for pagination parameters.

    Args:
        offset: Number of items to skip
        limit: Maximum number of items to return

    Returns:
        PaginationParams object with validated offset and limit
    """
    return PaginationParams(offset=offset, limit=limit)
