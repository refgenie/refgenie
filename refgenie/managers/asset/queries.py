"""
Shared SELECT builders for resolving ``genome/asset_group:asset_name``.

Every such lookup walks the same four tables (Asset, AssetName, AssetGroup,
Genome), so the join and its where-clause live here once; callers add their own
``.options(...)`` for the relationships they need eagerly loaded.

Keep this module dependent only on ``refgenie.db`` and SQLModel/SQLAlchemy.
Importing a manager here would make the manager/builder import cycle real.
"""

from collections.abc import Iterable

from sqlalchemy.orm import selectinload
from sqlalchemy.orm.interfaces import ORMOption
from sqlmodel import select
from sqlmodel.sql.expression import SelectOfScalar

from refgenie.db.tables import Asset, AssetGroup, AssetName, Genome, SeekKey
from refgenie.models import GenomeDigest


def asset_group_by_name_stmt(
    genome_digest: GenomeDigest, asset_group_name: str
) -> SelectOfScalar[AssetGroup]:
    """
    Select the asset group named ``genome_digest/asset_group_name``.

    Args:
        genome_digest: The digest of the genome.
        asset_group_name: The name of the asset group.

    Returns:
        The statement selecting the matching AssetGroup.
    """
    return (
        select(AssetGroup)
        .join(Genome)
        .where(AssetGroup.name == asset_group_name, Genome.digest == genome_digest)
    )


def asset_by_digest_stmt(digest: str) -> SelectOfScalar[Asset]:
    """
    Select the asset holding content ``digest``, with its seek keys and group loaded.

    Args:
        digest: The content digest of the asset.

    Returns:
        The statement selecting the matching Asset.
    """
    return (
        select(Asset)
        .where(Asset.digest == digest)
        .options(
            selectinload(Asset.seek_keys),
            selectinload(Asset.asset_group).selectinload(AssetGroup.genome),
        )
    )


def asset_by_name_stmt(
    genome_digest: GenomeDigest,
    asset_group_name: str,
    asset_name: str,
    *,
    options: Iterable[ORMOption] = (),
) -> SelectOfScalar[Asset]:
    """
    Select the asset named ``genome_digest/asset_group_name:asset_name``.

    Args:
        genome_digest: The digest of the genome.
        asset_group_name: The name of the asset group.
        asset_name: The name of the asset.
        options: Loader options (e.g. ``selectinload(...)``) to apply.

    Returns:
        The statement selecting the matching Asset.
    """
    return (
        select(Asset)
        .join(AssetName, AssetName.asset_digest == Asset.digest)
        .join(AssetGroup, AssetGroup.id == AssetName.asset_group_id)
        .join(Genome)
        .options(*options)
        .where(
            AssetName.name == asset_name,
            AssetGroup.name == asset_group_name,
            Genome.digest == genome_digest,
        )
    )


def seek_key_by_name_stmt(
    genome_digest: GenomeDigest,
    asset_group_name: str,
    asset_name: str,
    seek_key_name: str,
    *,
    options: Iterable[ORMOption] = (),
) -> SelectOfScalar[SeekKey]:
    """
    Select one named seek key of ``genome_digest/asset_group_name:asset_name``.

    Args:
        genome_digest: The digest of the genome.
        asset_group_name: The name of the asset group.
        asset_name: The name of the asset.
        seek_key_name: The name of the seek key.
        options: Loader options (e.g. ``selectinload(...)``) to apply.

    Returns:
        The statement selecting the matching SeekKey.
    """
    return (
        select(SeekKey)
        .join(Asset)
        .join(AssetName, AssetName.asset_digest == Asset.digest)
        .join(AssetGroup, AssetGroup.id == AssetName.asset_group_id)
        .join(Genome)
        .options(*options)
        .where(
            SeekKey.name == seek_key_name,
            AssetName.name == asset_name,
            AssetGroup.name == asset_group_name,
            Genome.digest == genome_digest,
        )
    )


def assets_stmt(genome_digests: Iterable[str] | None = None) -> SelectOfScalar[Asset]:
    """
    Select assets, with their group, genome, seek keys and names loaded.

    Args:
        genome_digests: Restrict to these genomes. None selects every asset.

    Returns:
        The statement selecting the assets. Callers may add further
        ``.where(...)`` clauses; ``AssetGroup`` and ``Genome`` are joined.
    """
    statement = (
        select(Asset)
        .join(AssetGroup)
        .join(Genome)
        .options(
            selectinload(Asset.asset_group).selectinload(AssetGroup.genome),
            selectinload(Asset.seek_keys),
            selectinload(Asset.asset_names),
        )
    )
    if genome_digests is not None:
        statement = statement.where(Genome.digest.in_(list(genome_digests)))
    return statement
