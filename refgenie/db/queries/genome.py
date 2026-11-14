"""Queries on the ``genome`` table."""

import json
from collections.abc import Collection, Mapping
from typing import Any

from sqlalchemy import JSON, Integer, String, bindparam, cast, false
from sqlalchemy import exists as sql_exists
from sqlalchemy.orm import selectinload
from sqlmodel import Session, func, or_, select

from refgenie.db.queries import species as species_queries
from refgenie.db.queries.common import one_or_raise
from refgenie.db.queries.paging import paginate
from refgenie.db.queries.search import (
    RankPatterns,
    like_escape,
    padded_words,
    rank_literal,
    search_conditions,
    validate_search_fields,
)
from refgenie.db.tables import AssetGroup, Genome
from refgenie.exceptions import MissingGenomeError
from refgenie.utils.search import (
    MatchRank,
    Search,
    SearchOperator,
    build_search_condition,
    search_words,
)

#: Fields a genome search may match, in the order the REST docs list them.
#: ``aliases`` is matched by the alias backend, not in SQL; the manager passes
#: the digests it found as ``alias_matches``. ``species`` is matched by
#: :func:`_species_condition`. ``synonyms`` (other names of the genome itself,
#: such as GRCh38) is a JSON list, matched on its text by :func:`_synonym_condition`.
SEARCH_FIELDS: tuple[str, ...] = (
    "digest",
    "description",
    "species",
    "synonyms",
    "assembly_source",
    "assembly_accession",
    "aliases",
)

#: Search fields that are not plain text columns of ``genome``.
_NON_TEXT_FIELDS = ("aliases", "synonyms", "species")

#: True for a genome with at least one asset group.
_HAS_ASSETS = sql_exists().where(AssetGroup.genome_digest == Genome.digest)

#: Asset groups of the genome in the outer query.
_ASSET_GROUP_COUNT = (
    select(func.count())
    .where(AssetGroup.genome_digest == Genome.digest)
    .correlate(Genome)
    .scalar_subquery()
)

#: Genomes with the most asset groups first, then by species, then by digest (unique).
_ORDER = (_ASSET_GROUP_COUNT.desc(), Genome.species_name.asc().nulls_last(), Genome.digest)


def get(session: Session, digest: str, *, with_aliases: bool = False) -> Genome:
    """
    The genome with this digest.

    Args:
        session: The session to execute in.
        digest: The genome digest.
        with_aliases: Also load the SQL ``aliases`` relationship.

    Raises:
        MissingGenomeError: If there is none.
    """
    statement = select(Genome).where(Genome.digest == digest)
    if with_aliases:
        statement = statement.options(selectinload(Genome.aliases))
    return one_or_raise(session, statement, MissingGenomeError(genome=digest))


def exists(session: Session, digest: str) -> bool:
    """Whether a genome with this digest exists."""
    return session.exec(select(Genome.digest).where(Genome.digest == digest)).first() is not None


def list_all(session: Session) -> list[Genome]:
    """Every genome, in listing order."""
    return list(session.exec(select(Genome).order_by(*_ORDER)).all())


def searches_aliases(search: Search | None) -> bool:
    """
    Whether ``search`` matches alias names, so the caller must find alias matches.

    Raises:
        ValueError: If the search names a field not in :data:`SEARCH_FIELDS`.
    """
    if search is None:
        return False
    validate_search_fields(search.fields, list(SEARCH_FIELDS))
    return bool(search.term.strip()) and "aliases" in (search.fields or SEARCH_FIELDS)


def _synonym_condition(term: str, operator: SearchOperator):
    """Match ``term`` against any one synonym, using the quotes around each in the JSON text."""
    text = cast(Genome.synonyms, String)
    # The JSON text escapes non-ASCII ("P\u00e8re"), so the term is escaped the
    # same way, and its backslashes doubled for LIKE.
    term = json.dumps(term)[1:-1].replace("\\", "\\\\")
    if operator == SearchOperator.EQUALS:
        return text.like(f'%"{term}"%', escape="\\")
    pattern = {
        SearchOperator.CONTAINS: f"%{term}%",
        SearchOperator.STARTS_WITH: f'%"{term}%',
        SearchOperator.ENDS_WITH: f'%{term}"%',
    }[operator]
    return text.ilike(pattern, escape="\\")


def _species_condition(term: str, operator: SearchOperator):
    """The genome's own species name, or any name of its taxon in the species list."""
    return or_(
        build_search_condition(Genome.species_name, term, operator),
        Genome.taxon_id.in_(species_queries.taxon_ids_named(term, operator)),
    )


#: Search fields whose matches are ranked by :class:`MatchRank`: the names of a
#: genome. A match in any other field (``description``) ranks ``OTHER``.
_RANKED_FIELDS = (
    "digest",
    "species",
    "synonyms",
    "assembly_source",
    "assembly_accession",
    "aliases",
)

#: Ranked fields that are plain text columns of ``genome``.
_RANKED_COLUMNS = ("digest", "assembly_source", "assembly_accession")


def _alias_ranks(alias_matches: Mapping[str, MatchRank], dialect: str):
    """``alias_matches`` as a subquery with columns ``digest`` and ``rank``.

    The digests travel as one JSON parameter, unpacked by ``json_each``, so
    the bind parameter count does not grow with the number of matches (a long
    ``IN`` list overflows the driver's limit on Postgres).
    """
    payload = bindparam(
        "alias_ranks", json.dumps({d: int(r) for d, r in alias_matches.items()}), type_=String
    )
    if dialect == "postgresql":
        pairs = func.json_each_text(cast(payload, JSON)).table_valued("key", "value")
        rank = cast(pairs.c.value, Integer)
    else:
        pairs = func.json_each(payload).table_valued("key", "value")
        rank = pairs.c.value
    return select(pairs.c.key.label("digest"), rank.label("rank")).subquery("alias_rank")


def _synonym_text(term: str) -> str:
    """``term`` as the JSON text of the ``synonyms`` column writes it, lowercased.

    The column is written by ``json.dumps``, which escapes non-ASCII ("P\\u00e8re"),
    so the term is escaped the same way. SQL ``lower`` then leaves the escapes'
    hex digits as they are, so non-ASCII letters compare case-sensitively, as
    in :func:`_synonym_condition`.
    """
    return json.dumps(term)[1:-1].lower()


def _ranked_search(
    term: str,
    fields,
    conditions: list,
    alias_ranks,
    dialect: str,
):
    """The genomes matching ``conditions``, each with its best :class:`MatchRank`.

    Returns:
        ``(matches, rank)``: a materialized CTE of matching digests with each
        ranked field padded once (see :func:`padded_words`), and the ORDER BY
        expression over its columns. None when the term has no words to rank
        by (only punctuation) or no ranked field is searched.
    """
    patterns = RankPatterns.for_term(term, "rank")
    ranked = [f for f in fields if f in _RANKED_FIELDS]
    if patterns is None or not ranked or not conditions:
        return None
    columns = [Genome.digest.label("digest")]
    statement = select(Genome)
    for field in ranked:
        if field in _RANKED_COLUMNS:
            columns.append(padded_words(getattr(Genome, field)).label(f"padded_{field}"))
        elif field == "species":
            columns.append(padded_words(Genome.species_name).label("species"))
            names = species_queries.taxon_name_ranks(patterns, search_words(term))
            statement = statement.outerjoin(names, names.c.taxon_id == Genome.taxon_id)
            columns.append(names.c.rank.label("species_list"))
        elif field == "synonyms":
            text = func.lower(cast(Genome.synonyms, String))
            columns.append(text.label("synonyms_text"))
            columns.append(padded_words(text).label("synonyms"))
        elif field == "aliases" and alias_ranks is not None:
            columns.append(alias_ranks.c.rank.label("aliases"))
    if alias_ranks is not None:
        statement = statement.outerjoin(alias_ranks, alias_ranks.c.digest == Genome.digest)
    matches = (
        statement.with_only_columns(*columns)
        .where(or_(*conditions))
        .cte("search_match")
        .prefix_with("MATERIALIZED")
    )

    other = rank_literal(MatchRank.OTHER)
    ranks = []
    for field in ranked:
        if field in _RANKED_COLUMNS:
            ranks.append(patterns.rank(matches.c[f"padded_{field}"]))
        elif field == "species":
            ranks.append(patterns.rank(matches.c.species))
            ranks.append(func.coalesce(matches.c.species_list, other))
        elif field == "synonyms":
            # Words of the JSON text would let "rat" equal one word of a
            # two-word synonym; an exact match is the whole quoted name.
            synonym = _synonym_text(term)
            synonym_patterns = RankPatterns.for_term(synonym, "rank_synonym")
            if synonym_patterns is not None:
                exact = matches.c.synonyms_text.like(
                    bindparam("rank_synonym_quoted", f'%"{like_escape(synonym)}"%'), escape="\\"
                )
                ranks.append(synonym_patterns.rank(matches.c.synonyms, exact=exact))
        elif field == "aliases" and alias_ranks is not None:
            ranks.append(func.coalesce(matches.c.aliases, other))
    if not ranks:
        return None
    best = func.least if dialect == "postgresql" else func.min
    return matches, ranks[0] if len(ranks) == 1 else best(*ranks)


def list_page(
    session: Session,
    *,
    digest: str | None = None,
    search: Search | None = None,
    alias_matches: Mapping[str, MatchRank] | None = None,
    has_assets: bool | None = None,
    taxon_id: int | None = None,
    digest_in: Collection[str] | None = None,
    served_stores: tuple[Collection[str], bool] | None = None,
    offset: int,
    limit: int,
) -> tuple[list[Genome], int]:
    """
    One page of genomes, most asset groups first, and the number that match.

    With a search term, better matches come first (:class:`MatchRank`: exact
    name, whole word of a name, word start of a name, anything else), and the
    usual order breaks ties.

    Args:
        session: The session to execute in.
        digest: Only this genome.
        search: Match ``search.term`` against the scalar search fields.
        alias_matches: Digests whose aliases matched ``search``, each with the
            best rank of its matching aliases; OR-ed into the search, since
            aliases may not live in SQL.
        has_assets: Only genomes with (True) or without (False) asset groups.
        taxon_id: Only genomes of this NCBI taxon.
        digest_in: Only genomes whose digest is in this set; an empty set
            matches nothing. ``None`` applies no digest-set filter.
        served_stores: ``(store names, include_storeless)``: only genomes
            whose ``store_name`` is one of these, or is unset when
            ``include_storeless``. ``None`` applies no store filter.
        offset: Rows to skip.
        limit: Most rows to return.

    Raises:
        ValueError: If the search names a field not in :data:`SEARCH_FIELDS`.
    """
    alias_matches = alias_matches or {}
    dialect = session.get_bind().dialect.name
    order_by: tuple[Any, ...] = _ORDER
    filters = []
    matches = None
    if digest is not None:
        filters.append(Genome.digest == digest)
    if search is not None:
        validate_search_fields(search.fields, list(SEARCH_FIELDS))
        if term := search.term.strip():
            fields = search.fields or SEARCH_FIELDS
            columns = [f for f in fields if f not in _NON_TEXT_FIELDS]
            conditions = search_conditions(Genome, term, search.operator, columns)
            alias_ranks = None
            if "aliases" in fields:
                if alias_matches:
                    alias_ranks = _alias_ranks(alias_matches, dialect)
                    conditions.append(Genome.digest.in_(select(alias_ranks.c.digest)))
                else:
                    conditions.append(false())
            if "synonyms" in fields:
                conditions.append(_synonym_condition(term, search.operator))
            if "species" in fields:
                conditions.append(_species_condition(term, search.operator))
            if conditions:
                search_filter = or_(*conditions)
                if ranked := _ranked_search(term, fields, conditions, alias_ranks, dialect):
                    matches, rank = ranked
                    order_by = (rank, *_ORDER)
                else:
                    filters.append(search_filter)
    if has_assets is not None:
        filters.append(_HAS_ASSETS if has_assets else ~_HAS_ASSETS)
    if taxon_id is not None:
        filters.append(Genome.taxon_id == taxon_id)
    if digest_in is not None:
        filters.append(Genome.digest.in_(list(digest_in)))
    if served_stores is not None:
        filters.append(served_condition(*served_stores))
    statement = select(Genome).where(*filters)
    if matches is None:
        return paginate(session, statement, order_by=order_by, offset=offset, limit=limit)
    # The page reads the ranks from the CTE; the count needs only the filter,
    # so it skips padding every name.
    return paginate(
        session,
        statement.join(matches, matches.c.digest == Genome.digest),
        count_statement=statement.where(search_filter),
        order_by=order_by,
        offset=offset,
        limit=limit,
    )


def served_condition(store_names: Collection[str], include_storeless: bool):
    """Genomes whose ``store_name`` is one of ``store_names`` (or unset, if asked)."""
    condition = Genome.store_name.in_(list(store_names))
    if include_storeless:
        condition = or_(condition, Genome.store_name.is_(None))
    return condition


def count_asset_groups(session: Session, digests: list[str]) -> dict[str, int]:
    """Asset groups per genome, for these digests, in one query. Absent means 0."""
    if not digests:
        return {}
    return dict(
        session.exec(
            select(AssetGroup.genome_digest, func.count(AssetGroup.id))
            .where(AssetGroup.genome_digest.in_(digests))
            .group_by(AssetGroup.genome_digest)
        ).all()
    )


def count_by_store(session: Session) -> dict[str, int]:
    """Genomes per owning store; genomes with no store are left out."""
    rows = session.exec(
        select(Genome.store_name, func.count(Genome.digest)).group_by(Genome.store_name)
    ).all()
    return {name: count for name, count in rows if name is not None}
