"""Parse ``?q=`` search query params from a request into SQL filters for the routers."""

from typing import Any

from fastapi import HTTPException, Query
from pydantic import BaseModel, Field
from sqlalchemy import Select, or_
from sqlmodel import SQLModel

from refgenie.utils.search import SearchOperator


class SearchParams(BaseModel):
    q: str | None = Field(None, description="General search query")
    search_fields: list[str] | None = Field(None, description="Fields to search in")
    operator: SearchOperator = Field(SearchOperator.CONTAINS, description="Search operator")


def _build_search_condition(column: Any, search_term: str, operator: SearchOperator) -> Any:
    """
    Build a search condition for a column based on the search operator.

    Args:
        column: SQLAlchemy column to search
        search_term: Search term to apply
        operator: Search operator to use

    Returns:
        SQLAlchemy condition
    """
    if operator == SearchOperator.EQUALS:
        return column == search_term
    elif operator == SearchOperator.CONTAINS:
        return column.ilike(f"%{search_term}%")
    elif operator == SearchOperator.STARTS_WITH:
        return column.ilike(f"{search_term}%")
    elif operator == SearchOperator.ENDS_WITH:
        return column.ilike(f"%{search_term}")


def _validate_search_fields(
    requested_fields: list[str] | None, searchable_columns: list[str]
) -> None:
    """
    Validate that all requested search fields are valid.

    Args:
        requested_fields: List of fields requested for searching
        searchable_columns: List of valid searchable column names

    Raises:
        HTTPException: If any requested fields are invalid
    """
    if not requested_fields:
        return

    if invalid_fields := [f for f in requested_fields if f not in searchable_columns]:
        raise HTTPException(
            status_code=422,
            detail=f"Invalid search fields: {', '.join(invalid_fields)}. "
            f"Valid fields are: {', '.join(searchable_columns)}",
        )


def _build_search_conditions(
    fields: list[str],
    search_term: str,
    operator: SearchOperator,
    model: type[SQLModel],
    relationship_fields: dict[str, str],
) -> list[Any]:
    """
    Build search conditions for the given fields.

    Args:
        fields: List of field names to search
        search_term: The search term to apply
        operator: The search operator to use
        model: The SQLModel class being queried
        relationship_fields: Dict mapping relationship field names to target fields

    Returns:
        List of SQLAlchemy conditions
    """
    search_conditions = []

    for field_name in fields:
        try:
            if field_name in relationship_fields:
                model_attr = getattr(model, field_name)
                related_model = model_attr.property.mapper.class_
                target_field = relationship_fields[field_name]

                if hasattr(related_model, target_field):
                    rel_column = getattr(related_model, target_field)
                    condition = _build_search_condition(rel_column, search_term, operator)
                    search_conditions.append(condition)
            else:
                model_attr = getattr(model, field_name)
                condition = _build_search_condition(model_attr, search_term, operator)
                search_conditions.append(condition)

        except (AttributeError, TypeError):
            # Skip fields that don't exist or can't be searched
            continue

    return search_conditions


def apply_search_to_query(
    query: Select[Any],
    search_params: SearchParams,
    searchable_columns: list[str],
    model: type[SQLModel],
    relationship_fields: dict[str, str] | None = None,
) -> Select[Any]:
    """
    Apply search parameters to a SQLModel query.

    Args:
        query: Base SQLModel query
        search_params: Search parameters
        searchable_columns: List of column names that can be searched
        model: The SQLModel class being queried
        relationship_fields: Dict mapping relationship field names to the target field
                            in the related model (e.g., {"aliases": "name"})

    Returns:
        Modified query with search filters applied
    """
    if not search_params.q:
        return query

    search_term = search_params.q.strip()
    if not search_term:
        return query

    _validate_search_fields(search_params.search_fields, searchable_columns)

    fields_to_search = search_params.search_fields or searchable_columns
    valid_fields = [f for f in fields_to_search if f in searchable_columns]

    if not valid_fields:
        return query

    relationship_fields = relationship_fields or {}
    search_conditions = _build_search_conditions(
        valid_fields,
        search_term,
        search_params.operator,
        model,
        relationship_fields,
    )

    if search_conditions:
        query = query.where(or_(*search_conditions))

    return query


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
