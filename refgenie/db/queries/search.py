"""The ``?q=`` text search, as a SQL condition over an entity's columns."""

from dataclasses import dataclass
from typing import Any

from sqlalchemy import Integer, String, bindparam, case, func, literal_column, or_
from sqlmodel import SQLModel

from refgenie.utils.search import (
    WORD_SEPARATORS,
    MatchRank,
    SearchOperator,
    build_search_condition,
    search_words,
)

#: Free-text fields: a ``contains`` search matches them only where a word
#: starts, so "rat" finds "rat genome" but not "curated" or "satellite".
#: Name fields (names, aliases, accessions) keep plain substring matching.
WORD_START_FIELDS = frozenset({"description"})


def _constant(text: str) -> Any:
    """``text`` written into the SQL, not bound: constants add no bind parameters."""
    return bindparam(None, text, type_=String, literal_execute=True)


def padded_words(column: Any) -> Any:
    """SQL for :func:`refgenie.utils.search.words`: ``column`` lowercased,
    separators as spaces, padded with a space each side.

    Plain ``lower``, ``replace`` and ``||``, so it runs on SQLite and Postgres
    alike (and faster on Postgres than one ``translate``). The separators are
    written into the SQL, so the chain adds no bind parameters.
    """
    expression = func.lower(column)
    for char in WORD_SEPARATORS:
        expression = func.replace(expression, _constant(char), _constant(" "))
    return _constant(" ") + expression + _constant(" ")


def like_escape(text: str) -> str:
    """``text`` with LIKE wildcards and the backslash escape char escaped."""
    return text.replace("\\", "\\\\").replace("%", "\\%").replace("_", "\\_")


def rank_condition(column: Any, term: str, rank: MatchRank) -> Any | None:
    """``column`` matches ``term`` at least as well as ``rank``; None for ``OTHER`` or a blank term.

    The SQL twin of :func:`refgenie.utils.search.match_rank`.
    """
    wanted = search_words(term)
    if not wanted or rank == MatchRank.OTHER:
        return None
    if rank == MatchRank.EXACT:
        return padded_words(column) == f" {wanted} "
    escaped = like_escape(wanted)
    pattern = f"% {escaped} %" if rank == MatchRank.WHOLE_WORD else f"% {escaped}%"
    return padded_words(column).like(pattern, escape="\\")


def word_start_condition(column: Any, term: str) -> Any:
    """``column`` has a word that starts with ``term``, case-insensitively.

    A term with no word characters (only punctuation) falls back to a substring match.
    """
    condition = rank_condition(column, term, MatchRank.WORD_START)
    return condition if condition is not None else column.ilike(f"%{term}%")


def rank_literal(rank: MatchRank) -> Any:
    """``rank`` written into the SQL as an integer, not bound."""
    return literal_column(str(int(rank)), Integer)


@dataclass(frozen=True)
class RankPatterns:
    """The :class:`MatchRank` tests of one search term, as bind parameters.

    Made once per term and shared by every column it ranks, so a named-parameter
    driver (psycopg) binds each pattern once however many columns read it.
    """

    exact: Any
    whole_word: Any
    word_start: Any

    @classmethod
    def for_term(cls, term: str, key: str) -> "RankPatterns | None":
        """The tests for ``term``, as parameters named after ``key``; None for a term with no words."""
        wanted = search_words(term)
        if not wanted:
            return None
        escaped = like_escape(wanted)
        return cls(
            exact=bindparam(f"{key}_exact", f" {wanted} "),
            whole_word=bindparam(f"{key}_whole", f"% {escaped} %"),
            word_start=bindparam(f"{key}_start", f"% {escaped}%"),
        )

    def rank(self, padded: Any, exact: Any = None) -> Any:
        """A CASE giving the :class:`MatchRank` of a :func:`padded_words` value.

        Args:
            padded: The padded value, best computed once (a column of a
                subquery), since the CASE reads it up to three times.
            exact: A test that replaces the ``EXACT`` one.
        """
        return case(
            (padded == self.exact if exact is None else exact, rank_literal(MatchRank.EXACT)),
            (padded.like(self.whole_word, escape="\\"), rank_literal(MatchRank.WHOLE_WORD)),
            (padded.like(self.word_start, escape="\\"), rank_literal(MatchRank.WORD_START)),
            else_=rank_literal(MatchRank.OTHER),
        )


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

    A ``contains`` search of a :data:`WORD_START_FIELDS` field matches word starts only.

    Returns:
        The conditions, to be OR-ed together. A field that does not exist on
        the model, or cannot be searched, is skipped.
    """
    conditions = []
    for field_name in fields:
        try:
            column = getattr(model, field_name)
            if operator == SearchOperator.CONTAINS and field_name in WORD_START_FIELDS:
                conditions.append(word_start_condition(column, term))
            else:
                conditions.append(build_search_condition(column, term, operator))
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
