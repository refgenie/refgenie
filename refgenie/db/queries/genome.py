"""Queries on the ``genome`` table."""

import json
from collections.abc import Iterable

from sqlalchemy import String, cast
from sqlalchemy import exists as sql_exists
from sqlalchemy.orm import selectinload
from sqlmodel import Session, func, or_, select

from refgenie.db.queries import species as species_queries
from refgenie.db.queries.common import one_or_raise
from refgenie.db.queries.paging import paginate
from refgenie.db.queries.search import search_conditions, validate_search_fields
from refgenie.db.tables import AssetGroup, Genome
from refgenie.exceptions import MissingGenomeError
from refgenie.utils.search import Search, SearchOperator, build_search_condition

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


def list_page(
    session: Session,
    *,
    digest: str | None = None,
    search: Search | None = None,
    alias_matches: Iterable[str] = (),
    has_assets: bool | None = None,
    taxon_id: int | None = None,
    offset: int,
    limit: int,
) -> tuple[list[Genome], int]:
    """
    One page of genomes, most asset groups first, and the number that match.

    Args:
        session: The session to execute in.
        digest: Only this genome.
        search: Match ``search.term`` against the scalar search fields.
        alias_matches: Digests whose aliases matched ``search``; OR-ed into
            the search, since aliases may not live in SQL.
        has_assets: Only genomes with (True) or without (False) asset groups.
        taxon_id: Only genomes of this NCBI taxon.
        offset: Rows to skip.
        limit: Most rows to return.

    Raises:
        ValueError: If the search names a field not in :data:`SEARCH_FIELDS`.
    """
    statement = select(Genome)
    if digest is not None:
        statement = statement.where(Genome.digest == digest)
    if search is not None:
        validate_search_fields(search.fields, list(SEARCH_FIELDS))
        if term := search.term.strip():
            fields = search.fields or SEARCH_FIELDS
            columns = [f for f in fields if f not in _NON_TEXT_FIELDS]
            conditions = search_conditions(Genome, term, search.operator, columns)
            if "aliases" in fields:
                conditions.append(Genome.digest.in_(list(alias_matches)))
            if "synonyms" in fields:
                conditions.append(_synonym_condition(term, search.operator))
            if "species" in fields:
                conditions.append(_species_condition(term, search.operator))
            statement = statement.where(or_(*conditions)) if conditions else statement
    if has_assets is not None:
        statement = statement.where(_HAS_ASSETS if has_assets else ~_HAS_ASSETS)
    if taxon_id is not None:
        statement = statement.where(Genome.taxon_id == taxon_id)
    return paginate(session, statement, order_by=_ORDER, offset=offset, limit=limit)


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
