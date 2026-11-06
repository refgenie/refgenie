"""Parse the ``?q=`` search query params of a list route into a ``Search``."""

from collections.abc import Iterator
from contextlib import contextmanager

from fastapi import HTTPException, Query
from pydantic import BaseModel, Field

from refgenie.utils.search import Search, SearchOperator


class SearchParams(BaseModel):
    q: str | None = Field(None, description="General search query")
    search_fields: list[str] | None = Field(None, description="Fields to search in")
    operator: SearchOperator = Field(SearchOperator.CONTAINS, description="Search operator")

    def to_search(self) -> Search | None:
        """The search the managers take, or None when ``q`` is blank."""
        if not (self.q or "").strip():
            return None
        return Search(term=self.q, fields=self.search_fields, operator=self.operator)


def get_search_params(
    q: str | None = Query(None, description="Search query"),
    search_fields: str | None = Query(None, description="Comma-separated list of fields to search"),
    operator: SearchOperator = Query(SearchOperator.CONTAINS, description="Search operator"),
) -> SearchParams:
    """FastAPI dependency for search parameters."""
    fields_list = None
    if search_fields:
        fields_list = [field.strip() for field in search_fields.split(",")]

    return SearchParams(
        q=q,
        search_fields=fields_list,
        operator=operator,
    )


@contextmanager
def invalid_search_is_422() -> Iterator[None]:
    """Turn the ``ValueError`` a manager raises for an unknown search field into a 422."""
    try:
        yield
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e)) from e
