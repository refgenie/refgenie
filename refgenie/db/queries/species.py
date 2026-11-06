"""Queries on the ``species`` and ``speciesname`` tables."""

from collections.abc import Iterable

from sqlmodel import Session, delete, select
from sqlalchemy.orm import selectinload

from refgenie.db.tables import Species, SpeciesName
from refgenie.utils.search import SearchOperator, build_search_condition


def replace_taxa(session: Session, species: list[Species]) -> None:
    """Rewrite these taxa: delete each one's row and names, then add ``species`` with their ``names``.

    Taxa not in ``species`` are not touched. ``species`` must have unique taxon ids.
    """
    ids = [s.taxon_id for s in species]
    if not ids:
        return
    session.exec(delete(SpeciesName).where(SpeciesName.taxon_id.in_(ids)))
    session.exec(delete(Species).where(Species.taxon_id.in_(ids)))
    session.add_all(species)


def list_by_taxon_ids(session: Session, taxon_ids: Iterable[int]) -> list[Species]:
    """The species with these taxon ids, names loaded, by taxon id. Ids with no row are left out."""
    ids = list(taxon_ids)
    if not ids:
        return []
    statement = (
        select(Species)
        .where(Species.taxon_id.in_(ids))
        .options(selectinload(Species.names))
        .order_by(Species.taxon_id)
    )
    return list(session.exec(statement).all())


def taxon_ids_named(term: str, operator: SearchOperator):
    """The SELECT of taxon ids with any name (scientific, common or other) matching ``term``.

    A statement, not rows: the genome ``species`` search field uses it as a
    subquery, ``Genome.taxon_id.in_(taxon_ids_named(...))``. It is the one copy
    of the species matching rules.
    """
    return (
        select(SpeciesName.taxon_id)
        .where(build_search_condition(SpeciesName.name, term, operator))
        .distinct()
    )


def matching_taxon_ids(session: Session, term: str, operator: SearchOperator) -> list[int]:
    """Taxon ids with any name matching ``term`` under ``operator``, ascending. A blank term matches nothing."""
    term = term.strip()
    if not term:
        return []
    return list(session.exec(taxon_ids_named(term, operator).order_by(SpeciesName.taxon_id)).all())
