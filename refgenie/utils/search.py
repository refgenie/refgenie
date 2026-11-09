"""Search operators and the match rules shared by the SQL search and the in-memory alias search."""

from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum, IntEnum
from typing import Any


class SearchOperator(str, Enum):
    EQUALS = "eq"
    CONTAINS = "contains"
    STARTS_WITH = "starts_with"
    ENDS_WITH = "ends_with"


@dataclass(frozen=True)
class Search:
    """A text search over an entity's searchable fields, as managers take it.

    ``fields`` None searches every searchable field. A blank ``term`` filters
    nothing.
    """

    term: str
    fields: list[str] | None = None
    operator: SearchOperator = SearchOperator.CONTAINS


class MatchRank(IntEnum):
    """How well a name matches a search term; lower ranks list first.

    Names compare case-insensitively, with punctuation (:data:`WORD_SEPARATORS`)
    read as spaces, so "Norway rat" holds the whole word "rat" and
    "mRatBN7.2" holds "rat" only inside a word.
    """

    #: The whole name is the term ("human" for the common name "human").
    EXACT = 0
    #: The term is a whole word of the name ("rat" in "Norway rat").
    WHOLE_WORD = 1
    #: A word of the name starts with the term ("rat" in "Rattus norvegicus").
    WORD_START = 2
    #: Any other match: inside a word, or in a field that is not ranked.
    OTHER = 3


#: Characters read as spaces when a search splits text into words.
WORD_SEPARATORS = "-_()[]{}/\\,.:;\"'\t\n\r"


def words(text: str) -> str:
    """``text`` lowercased, separators as spaces, padded with one space each side.

    The Python twin of the SQL ``padded_words`` in ``refgenie.db.queries.search``:
    a whole word ``w`` is then the substring `` w ``, a word start `` w``.
    """
    lowered = text.lower()
    for char in WORD_SEPARATORS:
        lowered = lowered.replace(char, " ")
    return f" {lowered} "


def search_words(term: str) -> str:
    """The term as :func:`words` reads it, without the padding and with single spaces."""
    return " ".join(words(term).split())


def match_rank(value: str, term: str) -> MatchRank:
    """The :class:`MatchRank` of a name that matched ``term``."""
    padded, wanted = words(value), search_words(term)
    if not wanted:
        return MatchRank.OTHER
    if padded.strip() == wanted:
        return MatchRank.EXACT
    if f" {wanted} " in padded:
        return MatchRank.WHOLE_WORD
    if f" {wanted}" in padded:
        return MatchRank.WORD_START
    return MatchRank.OTHER


def build_search_condition(column: Any, search_term: str, operator: SearchOperator) -> Any:
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
    Python-side equivalent of ``build_search_condition`` for in-memory matching.

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
) -> Mapping[str, MatchRank]:
    """
    Genome digests whose alias names match a search term, each with its best match rank.

    Reads aliases from the mode-selected alias manager (SQL-backed in server mode,
    store-backed in local mode) so alias search works in both modes without a SQL join.

    Args:
        alias_manager: The alias manager (``refgenie.alias``); records expose
            ``.name`` and ``.genome_digest``.
        search_term: The search term to apply to alias names.
        operator: The search operator to use.

    Returns:
        Matching genome digest -> the best :func:`match_rank` of its matching aliases.
    """
    ranks: dict[str, MatchRank] = {}
    for record in alias_manager.list_all():
        if name_matches_operator(record.name, search_term, operator):
            rank = match_rank(record.name, search_term)
            ranks[record.genome_digest] = min(rank, ranks.get(record.genome_digest, rank))
    return ranks
