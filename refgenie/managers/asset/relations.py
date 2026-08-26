"""
AssetRelations - manages asset parent-child relationships and provenance.

Internal helper class used by builder, puller, and remove operations.
Not user-facing - users interact with AssetManager instead.
"""


from sqlalchemy.engine import Engine
from sqlalchemy.orm import selectinload
from sqlmodel import select

from refgenie.db.tables import Asset, AssetGroup
from refgenie.logger import logger
from refgenie.managers.asset.queries import asset_by_name_stmt
from refgenie.managers.base import ResourceManager


class AssetRelations(ResourceManager):
    """
    Manages asset parent-child relationships and provenance.

    This is an internal helper class - not exposed to users.
    Used by builder, puller, and remove operations to manage
    asset dependency relationships.
    """

    def __init__(self, database_engine: Engine):
        """
        Initialize the AssetRelations manager.

        Args:
            database_engine: The database engine.
        """
        super().__init__(database_engine)

    def get_parents(
        self, genome_digest: str, asset_group_name: str, asset_name: str
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
            asset = (
                session.exec(asset_by_name_stmt(genome_digest, asset_group_name, asset_name))
                .unique()
                .one()
            )
            return asset.parents

    def get_children(
        self,
        asset_group_name: str,
        asset_name: str,
        genome_digest: str,
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
            asset = (
                session.exec(
                    asset_by_name_stmt(
                        genome_digest,
                        asset_group_name,
                        asset_name,
                        options=[selectinload(Asset.children)],
                    )
                )
                .unique()
                .one()
            )
            # Now eagerly load the asset_group and genome for each child
            children_with_relationships = []
            for child in asset.children:
                child_with_rel = (
                    session.exec(
                        select(Asset)
                        .options(selectinload(Asset.asset_group).selectinload(AssetGroup.genome))
                        .where(Asset.digest == child.digest)
                    )
                    .unique()
                    .one()
                )
                children_with_relationships.append(child_with_rel)

            return children_with_relationships

    def get_size(
        self,
        asset_group_name: str,
        asset_name: str,
        genome_digest: str,
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
            asset = (
                session.exec(
                    asset_by_name_stmt(
                        genome_digest,
                        asset_group_name,
                        asset_name,
                        options=[selectinload(Asset.seek_keys)],
                    )
                )
                .unique()
                .one()
            )
            return sum(seek_key.size for seek_key in asset.seek_keys if seek_key.size)

    def set_parents(
        self,
        genome_digest: str,
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
            asset = (
                session.exec(asset_by_name_stmt(genome_digest, asset_group_name, asset_name))
                .unique()
                .one()
            )
            parent_assets = (
                session.exec(select(Asset).where(Asset.digest.in_(parent_asset_digests)))
                .unique()
                .all()
            )
            asset.parents = list(parent_assets)
            session.commit()
        logger.info(f"Updated parents of '{genome_digest}/{asset_group_name}:{asset_name}'")

    def set_children(
        self,
        genome_digest: str,
        asset_group_name: str,
        asset_name: str,
        child_asset_digests: list[str],
    ):
        """
        Set the children assets of an asset.

        Args:
            genome_digest: The digest of the genome.
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            child_asset_digests: The digests of the child assets.
        """
        with self._database_session as session:
            asset = (
                session.exec(asset_by_name_stmt(genome_digest, asset_group_name, asset_name))
                .unique()
                .one()
            )
            child_assets = (
                session.exec(select(Asset).where(Asset.digest.in_(child_asset_digests)))
                .unique()
                .all()
            )
            asset.children = list(child_assets)
            session.commit()
        logger.info(f"Updated children of asset '{asset_name}'")

    def parent_digests_exist(self, parent_assets: dict, asset_exists_callback) -> bool:
        """
        Check if all parent assets exist in the database.

        Args:
            parent_assets: A list of parent assets, each with 'name' and 'digest'.
            asset_exists_callback: Callback to check if asset exists by digest.

        Returns:
            bool: True if all parent assets exist, False otherwise.
        """
        digests = [(asset["name"], asset["digest"]) for asset in parent_assets]
        for name, digest in digests:
            if not asset_exists_callback(digest=digest):
                logger.error(
                    f"Parent asset with digest {digest} does not exist. Server asset name: {name}"
                )
                return False
        return True
