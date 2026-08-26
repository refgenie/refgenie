from enum import Enum
from typing import Any

from fastapi import HTTPException, Query
from pydantic import BaseModel, Field
from sqlalchemy import Select, or_
from sqlmodel import SQLModel


class SearchOperator(str, Enum):
    EQUALS = "eq"
    CONTAINS = "contains"
    STARTS_WITH = "starts_with"
    ENDS_WITH = "ends_with"


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


def name_matches_operator(value: str, search_term: str, operator: SearchOperator) -> bool:
    """
    Python-side equivalent of ``_build_search_condition`` for in-memory matching.

    Mirrors the SQL operator semantics: EQUALS is case-sensitive; CONTAINS,
    STARTS_WITH and ENDS_WITH are case-insensitive (matching SQL ``ilike``).

    Args:
        value: The string value to test (e.g. an alias name).
        search_term: The search term to apply.
        operator: The search operator to use.

    Returns:
        True if ``value`` matches ``search_term`` under ``operator``.
    """
    if operator == SearchOperator.EQUALS:
        return value == search_term
    lowered_value = value.lower()
    lowered_term = search_term.lower()
    if operator == SearchOperator.CONTAINS:
        return lowered_term in lowered_value
    elif operator == SearchOperator.STARTS_WITH:
        return lowered_value.startswith(lowered_term)
    elif operator == SearchOperator.ENDS_WITH:
        return lowered_value.endswith(lowered_term)
    return False


def alias_digests_matching(
    alias_manager: Any, search_term: str, operator: SearchOperator
) -> list[str]:
    """
    Return genome digests whose alias names match a search term via the alias manager.

    Reads aliases from the mode-selected alias manager (SQL-backed in server mode,
    store-backed in local mode) so alias search works in both modes without a SQL join.

    Args:
        alias_manager: The alias manager (``refgenie.alias``); records expose
            ``.name`` and ``.genome_digest``.
        search_term: The search term to apply to alias names.
        operator: The search operator to use.

    Returns:
        List of matching genome digests (may contain duplicates).
    """
    return [
        record.genome_digest
        for record in alias_manager.list_all()
        if name_matches_operator(record.name, search_term, operator)
    ]


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
    search_fields: str | None = Query(
        None, description="Comma-separated list of fields to search"
    ),
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
