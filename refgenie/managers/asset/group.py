"""
AssetGroupManager - asset groups and the default asset of each group.

Reached as ``rgc.asset.group``. Needs only the database and the genome manager; it
never calls back into asset code, so content, seek keys, the alias tree, the
builder and the puller can all depend on it.
"""

from collections.abc import Iterable
from typing import TYPE_CHECKING

from sqlalchemy import update as sa_update
from sqlalchemy.engine import Engine
from sqlalchemy.orm import selectinload
from sqlmodel import Session, select

from refgenie.db.tables import AssetGroup, AssetName, Genome
from refgenie.exceptions import MissingAssetError, MissingAssetGroupError
from refgenie.logger import logger
from refgenie.managers.asset.queries import asset_group_by_name_stmt
from refgenie.managers.base import ResourceManager
from refgenie.managers.queries import one_or_raise
from refgenie.models import GenomeDigest
from refgenie.plugins.events import NULL_EVENTS, EventSink, update_scope
from refgenie.plugins.hooks import Change

if TYPE_CHECKING:
    from refgenie.managers.genome import GenomeManager


def asset_group_in_session(
    session: Session,
    asset_group_name: str,
    *,
    genome_digest: GenomeDigest,
) -> AssetGroup:
    """
    ``AssetGroupManager.get``, asked of an already-open session.

    For callers that hold a session: ``_database_session`` must not be
    re-entered from inside an open block.

    Raises:
        MissingAssetGroupError: If the asset group does not exist.
    """
    return one_or_raise(
        session,
        asset_group_by_name_stmt(genome_digest, asset_group_name),
        MissingAssetGroupError(genome=genome_digest, asset_group=asset_group_name),
        unique=True,
    )


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
            return asset_group_in_session(session, asset_group_name, genome_digest=genome_digest)

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
            result = session.exec(asset_group_by_name_stmt(genome_digest, asset_group_name))
            return bool(result.first())

    def list_all(self, genome_digests: list[GenomeDigest] | None = None) -> Iterable[AssetGroup]:
        """
        List all asset groups.

        Args:
            genome_digests: Genome digests to filter by.

        Returns:
            Iterable[AssetGroup]: A list of all asset groups.
        """
        statement = (
            select(AssetGroup)
            .join(Genome)
            .options(selectinload(AssetGroup.genome))
            .options(selectinload(AssetGroup.assets))
        )
        if genome_digests is not None:
            statement = statement.where(Genome.digest.in_(genome_digests))

        with self._database_session as session:
            result = session.exec(statement)
            return result.unique().all()

    @update_scope
    def remove_rows(self, asset_group_name: str, genome_digest: GenomeDigest):
        """
        Delete an asset group's rows, and its genome's if it was the last group.

        Rows only: no disk cleanup. The group's ``alias/`` and ``builds/`` trees
        stay on disk. Callers normally want ``rgc.asset.remove``, which removes
        an asset and, with its last asset, the group and its trees
        (``AssetManager.remove_by_digest`` is what calls this).

        Args:
            asset_group_name: The name of the asset group.
            genome_digest: The digest of the genome.
        """
        with self._database_session as session:
            asset_group = asset_group_in_session(
                session, asset_group_name, genome_digest=genome_digest
            )
            session.delete(asset_group)
            session.commit()
            genome_is_empty = (
                session.exec(
                    select(AssetGroup).join(Genome).where(Genome.digest == genome_digest)
                ).first()
                is None
            )
        # The genome removal is a write in another manager, so it runs after this
        # session closes rather than opening a second one inside it.
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
            # Look up the exact AssetName row for the name the caller passed,
            # scoped to this group. A name from another group is simply not found
            # here -- that is the same-group guard, by construction. The name the
            # caller chose is what gets flagged, so a non-canonical name is not
            # silently canonicalized.
            target = one_or_raise(
                session,
                select(AssetName).where(
                    AssetName.asset_group_id == asset_group.id,
                    AssetName.name == asset_name,
                ),
                MissingAssetError(
                    genome=genome_digest, asset_group=asset_group_name, asset=asset_name
                ),
            )
            previous = session.exec(
                select(AssetName.name).where(
                    AssetName.asset_group_id == asset_group.id,
                    AssetName.is_default == True,  # noqa: E712 - SQL boolean, not Python identity
                )
            ).one_or_none()
            # Clear the current default first, then set the new one. The partial
            # unique index forbids two defaults; a single transaction makes the
            # intermediate zero-default state invisible.
            session.exec(
                sa_update(AssetName)
                .where(
                    AssetName.asset_group_id == asset_group.id,
                    AssetName.is_default == True,  # noqa: E712 - SQL boolean, not Python identity
                )
                .values(is_default=False)
            )
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
            row = session.exec(
                select(AssetName).where(
                    AssetName.asset_group_id == asset_group.id,
                    AssetName.is_default == True,  # noqa: E712 - SQL boolean, not Python identity
                )
            ).one_or_none()
            return row.name if row else None
