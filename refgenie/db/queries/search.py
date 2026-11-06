"""The ``?q=`` text search, as a SQL condition over an entity's columns."""

from typing import Any

from sqlalchemy import or_
from sqlmodel import SQLModel

from refgenie.utils.search import SearchOperator, build_search_condition


def validate_search_fields(requested: list[str] | None, searchable: list[str]) -> None:
    """
    Raise ``ValueError`` if any requested field is not searchable.

    Args:
        requested: The fields the caller asked to search, or None for all.
        searchable: The entity's searchable field names.
    """
    if requested and (invalid := [f for f in requested if f not in searchable]):
        raise ValueError(
            f"Invalid search fields: {', '.join(invalid)}. "
            f"Valid fields are: {', '.join(searchable)}"
        )


def search_conditions(
    model: type[SQLModel], term: str, operator: SearchOperator, fields: list[str]
) -> list[Any]:
    """
    One condition per field: ``term`` matched against that column of ``model``.

    Returns:
        The conditions, to be OR-ed together. A field that does not exist on
        the model, or cannot be searched, is skipped.
    """
    conditions = []
    for field_name in fields:
        try:
            conditions.append(build_search_condition(getattr(model, field_name), term, operator))
        except (AttributeError, TypeError):
            continue
    return conditions


def search_condition(
    model: type[SQLModel],
    term: str | None,
    operator: SearchOperator,
    fields: list[str] | None,
    searchable: list[str],
) -> Any | None:
    """
    The OR of every searched field matching ``term``, or None for no filter.

    Args:
        model: The table class being queried.
        term: The search term; blank means no filter.
        operator: How the term matches.
        fields: The fields to search, or None for every searchable one.
        searchable: The entity's searchable field names.

    Raises:
        ValueError: If ``fields`` names a field not in ``searchable``.
    """
    validate_search_fields(fields, searchable)
    term = (term or "").strip()
    if not term:
        return None
    wanted = [f for f in (fields or searchable) if f in searchable]
    conditions = search_conditions(model, term, operator, wanted)
    return or_(*conditions) if conditions else None
