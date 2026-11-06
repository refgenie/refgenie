"""
AssetLinkManager - parent/child links between assets.

Reached as ``rgc.asset.links``. The builder records the inputs a build used as
the new asset's parents; the puller copies the parents a server reports.
"""

from refgenie.db.engine import Engine
from refgenie.db.queries import asset as queries
from refgenie.db.tables import Asset
from refgenie.logger import logger
from refgenie.managers.base import ResourceManager
from refgenie.models import GenomeDigest


class AssetLinkManager(ResourceManager):
    """Parent/child links between assets."""

    def __init__(self, database_engine: Engine):
        """
        Initialize the AssetLinkManager.

        Args:
            database_engine: The database engine.
        """
        super().__init__(database_engine)

    def parents(
        self, genome_digest: GenomeDigest, asset_group_name: str, asset_name: str
    ) -> list[Asset]:
        """
        Get the parent assets of an asset.

        Args:
            genome_digest: The digest of the genome.
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.

        Returns:
            list[Asset]: The parent assets of the asset.
        """
        with self._database_session as session:
            return queries.get_by_name(
                session, genome_digest, asset_group_name, asset_name, with_links=True
            ).parents

    def children(
        self,
        asset_group_name: str,
        asset_name: str,
        genome_digest: GenomeDigest,
    ) -> list[Asset]:
        """
        Get the children assets of an asset.

        Args:
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            genome_digest: The digest of the genome.

        Returns:
            list[Asset]: The children assets of the asset with eager-loaded relationships.
        """
        with self._database_session as session:
            return queries.get_by_name(
                session, genome_digest, asset_group_name, asset_name, with_links=True
            ).children

    def size(
        self,
        asset_group_name: str,
        asset_name: str,
        genome_digest: GenomeDigest,
    ) -> int:
        """
        Get the size of an asset in bytes.

        Args:
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            genome_digest: The digest of the genome.

        Returns:
            int: The size of the asset (sum of all seek key sizes).
        """
        with self._database_session as session:
            asset = queries.get_by_name(session, genome_digest, asset_group_name, asset_name)
            return sum(seek_key.size for seek_key in asset.seek_keys if seek_key.size)

    def set_parents(
        self,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        asset_name: str,
        parent_asset_digests: list[str],
    ):
        """
        Update the parent assets of an asset.

        Args:
            genome_digest: The digest of the genome.
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            parent_asset_digests: The digests of the parent assets.
        """
        with self._database_session as session:
            asset = queries.get_by_name(session, genome_digest, asset_group_name, asset_name)
            parent_assets = queries.list_by_digests(session, parent_asset_digests)
            asset.parents = list(parent_assets)
            session.commit()
        logger.info(f"Updated parents of '{genome_digest}/{asset_group_name}:{asset_name}'")
