"""Catalog-wide counts: the summary endpoints and the ping probe."""

from sqlmodel import Session, func, select

from refgenie.db.tables import Asset, AssetClass, AssetGroup, Genome


def count_all(session: Session) -> dict[str, int]:
    """How many genomes, asset groups and assets there are."""
    return {
        "genomes": session.exec(select(func.count()).select_from(Genome)).one(),
        "asset_groups": session.exec(select(func.count()).select_from(AssetGroup)).one(),
        "assets": session.exec(select(func.count()).select_from(Asset)).one(),
    }


def by_species(session: Session) -> list[tuple[str | None, int, int, int]]:
    """
    ``(species, genomes, asset classes, assets)`` for every species with assets.

    Only genomes with at least one asset count, since the joins are inner.
    Rows come by species name, the unnamed species last.
    """
    return list(
        session.exec(
            select(
                Genome.species_name,
                func.count(func.distinct(Genome.digest)),
                func.count(func.distinct(AssetClass.id)),
                func.count(Asset.digest),
            )
            .select_from(Genome)
            .join(AssetGroup, Genome.digest == AssetGroup.genome_digest)
            .join(AssetClass, AssetGroup.asset_class_id == AssetClass.id)
            .join(Asset, AssetGroup.id == Asset.asset_group_id)
            .group_by(Genome.species_name)
            .order_by(Genome.species_name.asc().nulls_last())
        ).all()
    )
