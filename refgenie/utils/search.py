"""Search operators and the match rules shared by the SQL search and the in-memory alias search."""

from enum import Enum
from typing import Any


class SearchOperator(str, Enum):
    EQUALS = "eq"
    CONTAINS = "contains"
    STARTS_WITH = "starts_with"
    ENDS_WITH = "ends_with"


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
