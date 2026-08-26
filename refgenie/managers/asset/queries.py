"""
Shared SELECT builders for resolving ``genome/asset_group:asset_name``.

Every such lookup walks the same four tables (Asset, AssetName, AssetGroup,
Genome), so the join and its where-clause live here once; callers add their own
``.options(...)`` for the relationships they need eagerly loaded.

Keep this module dependent only on ``refgenie.db`` and SQLModel/SQLAlchemy.
Importing a manager here would make the manager/builder import cycle real.
"""

from collections.abc import Iterable

from sqlalchemy.orm.interfaces import ORMOption
from sqlmodel import select
from sqlmodel.sql.expression import SelectOfScalar

from refgenie.db.tables import Asset, AssetGroup, AssetName, Genome, SeekKey


def asset_by_name_stmt(
    genome_digest: str,
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
    genome_digest: str,
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
