"""
AssetGroupManager - asset groups and the default asset of each group.

Reached as ``rgc.asset.group``. Needs only the database and the genome manager; it
never calls back into asset code, so content, seek keys, the alias tree, the
builder and the puller can all depend on it.
"""

from collections.abc import Iterable
from typing import TYPE_CHECKING

from refgenie.db.engine import Engine
from refgenie.db.queries import asset_group as queries
from refgenie.db.tables import AssetGroup
from refgenie.exceptions import MissingAssetError
from refgenie.logger import logger
from refgenie.managers.base import ResourceManager
from refgenie.models import GenomeDigest, Page
from refgenie.plugins.events import NULL_EVENTS, EventSink, update_scope
from refgenie.plugins.hooks import Change
from refgenie.utils.search import Search

if TYPE_CHECKING:
    from refgenie.managers.genome import GenomeManager


class AssetGroupManager(ResourceManager):
    """Manager for asset groups and their default asset."""

    def __init__(
        self,
        database_engine: Engine,
        genome_manager: "GenomeManager",
        events: EventSink | None = None,
    ):
        """
        Initialize the AssetGroupManager.

        Args:
            database_engine: The database engine.
            genome_manager: The GenomeManager, for removing a genome whose last
                group is removed.
            events: Where default changes are recorded for plugins.
        """
        super().__init__(database_engine)
        self._genome_manager = genome_manager
        self._events = events or NULL_EVENTS

    def get(
        self,
        asset_group_name: str,
        *,
        genome_digest: GenomeDigest,
    ) -> AssetGroup:
        """
        Get an asset group by its name.

        Args:
            asset_group_name: The name of the asset group.
            genome_digest: The genome digest.

        Returns:
            AssetGroup: The asset group.
        """
        with self._database_session as session:
            return queries.get(session, genome_digest, asset_group_name)

    def get_by_id(self, id: int) -> AssetGroup | None:
        """The asset group with this id, or None if there is none."""
        with self._database_session as session:
            return queries.get_by_id(session, id)

    def exists(
        self,
        asset_group_name: str,
        *,
        genome_digest: GenomeDigest,
    ) -> bool:
        """
        Check if an asset group exists.

        Args:
            asset_group_name: The name of the asset group.
            genome_digest: The genome digest.

        Returns:
            bool: Whether the asset group exists. False for an unknown genome.
        """
        with self._database_session as session:
            return queries.exists(session, genome_digest, asset_group_name)

    def list_all(self, genome_digests: list[GenomeDigest] | None = None) -> Iterable[AssetGroup]:
        """
        List all asset groups.

        Args:
            genome_digests: Genome digests to filter by.

        Returns:
            Iterable[AssetGroup]: A list of all asset groups.
        """
        with self._database_session as session:
            return queries.list_all(session, genome_digests)

    def page(
        self,
        *,
        genome_digest: GenomeDigest | None = None,
        asset_class: str | None = None,
        name: str | None = None,
        id: int | None = None,
        search: Search | None = None,
        offset: int,
        limit: int,
    ) -> Page[AssetGroup]:
        """
        One page of asset groups, by genome, name and id.

        Raises:
            ValueError: If the search names a field other than ``name``.
        """
        with self._database_session as session:
            rows, total = queries.list_page(
                session,
                genome_digest=genome_digest,
                asset_class=asset_class,
                name=name,
                id=id,
                search=search,
                offset=offset,
                limit=limit,
            )
        return Page(items=rows, offset=offset, limit=limit, total=total)

    @update_scope
    def remove_rows(self, asset_group_name: str, genome_digest: GenomeDigest):
        """
        Delete an asset group's rows, and its genome's if it was the last group.

        Rows only: the ``alias/`` and ``builds/`` trees stay on disk. Callers
        normally want ``rgc.asset.remove``, which also cleans those up.

        Args:
            asset_group_name: The name of the asset group.
            genome_digest: The digest of the genome.
        """
        with self._database_session as session:
            asset_group = queries.get(session, genome_digest, asset_group_name)
            session.delete(asset_group)
            session.commit()
            genome_is_empty = not queries.exists_for_genome(session, genome_digest)
        # After the session closes: the genome removal opens its own.
        if genome_is_empty:
            logger.info(
                f"No assets groups left for the genome '{genome_digest}'. Removing the genome"
            )
            self._genome_manager.remove(genome_digest)
        logger.info(f"Removed asset group and all assets '{genome_digest}/{asset_group_name}'")

    @update_scope
    def set_default(
        self,
        asset_group_name: str,
        asset_name: str,
        *,
        genome_digest: GenomeDigest,
    ):
        """
        Set the default asset for an asset group.

        Args:
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            genome_digest: The genome digest.
        """
        asset_group = self.get(asset_group_name, genome_digest=genome_digest)
        with self._database_session as session:
            # Scoped to this group, so a name from another group is not found.
            target = queries.find_name(session, asset_group.id, asset_name)
            if target is None:
                raise MissingAssetError(
                    genome=genome_digest, asset_group=asset_group_name, asset=asset_name
                )
            previous = queries.get_default_name(session, asset_group.id)
            # Clear first: a partial unique index forbids two defaults.
            queries.clear_default(session, asset_group.id)
            target.is_default = True
            session.add(target)
            session.commit()
            asset_digest = target.asset_digest
        logger.info(f"Set default asset: '{genome_digest}/{asset_group_name}:{asset_name}'")
        if previous != asset_name:
            self._events.record(
                Change(
                    action="default_changed",
                    genome=genome_digest,
                    asset_group=asset_group_name,
                    asset=asset_name,
                    previous=previous,
                    digest=asset_digest,
                )
            )

    def get_default(
        self,
        asset_group_name: str,
        *,
        genome_digest: GenomeDigest,
    ) -> str | None:
        """
        Get the default asset for an asset group.

        Args:
            asset_group_name: The name of the asset group.
            genome_digest: The genome digest.

        Returns:
            str | None: The name of the default asset, or None if no default is set.

        Raises:
            MissingAssetGroupError: If the asset group does not exist.
        """
        asset_group = self.get(asset_group_name, genome_digest=genome_digest)
        with self._database_session as session:
            return queries.get_default_name(session, asset_group.id)
