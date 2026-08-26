"""SQLAlchemy event handlers.

The destructive handlers here do not destroy anything. They resolve the paths a
delete makes obsolete -- which can only be done while the row is still present
-- and hand them to ``refgenie.db.cleanup``, which runs them after the
transaction commits. Session-level listeners drain that queue on ``after_commit``
and discard it on rollback, so a rolled-back delete leaves the files intact
alongside the rows it restored.
"""

from pathlib import Path

from sqlalchemy import Connection, event, func
from sqlalchemy.orm import Session
from sqlalchemy.sql import select
from refgenie.db.cleanup import (
    discard_pending_cleanups,
    enqueue_callable,
    enqueue_path_removal,
    run_pending_cleanups,
)
from refgenie.db.tables import (
    Alias,
    StagedAsset,
    Asset,
    AssetGroup,
    Configuration,
    SeekKey,
    SeekKeyType,
)
from refgenie.logger import logger
from refgenie.utils.symlinks import alias_owned_paths


def _genome_folder(connection: Connection) -> str:
    """Resolve the configured genome folder from the newest Configuration row."""
    return connection.execute(
        select(Configuration.genome_folder).where(
            Configuration.id == select(func.max(Configuration.id)).scalar_subquery()
        )
    ).scalar_one()

_events_registered = False


def register_events():
    """Register all SQLAlchemy event handlers. Safe to call multiple times."""
    global _events_registered
    if _events_registered:
        return

    event.listen(Asset, "before_delete", _asset_before_delete_handler)
    event.listen(Alias, "before_delete", _alias_before_delete_handler)
    event.listen(StagedAsset, "before_delete", _staged_asset_before_delete_handler)
    event.listen(SeekKey, "before_insert", _seek_key_before_insert_handler)

    # Session-level: the deferred cleanup queue. Listening on the base Session
    # class covers every subclass, including SQLModel's.
    event.listen(Session, "after_commit", _session_after_commit_handler)
    event.listen(Session, "after_rollback", _session_after_rollback_handler)
    event.listen(Session, "after_soft_rollback", _session_after_soft_rollback_handler)

    _events_registered = True


def _session_after_commit_handler(session: Session):
    """Run the filesystem work the just-committed transaction authorized.

    This must not emit SQL: the transaction is over. Filesystem work only.
    """
    run_pending_cleanups(session)


def _session_after_rollback_handler(session: Session):
    discard_pending_cleanups(session)


def _session_after_soft_rollback_handler(session: Session, _previous_transaction):
    discard_pending_cleanups(session)


def _asset_before_delete_handler(_, connection: Connection, target: Asset):
    """
    Queue the asset's content directory for removal after the commit.

    The content directory is digest-addressed and shared by every name the
    content carries; removing the asset removes the content. Alias-tree
    directories (one per name) are cleared by the AssetManager.remove path,
    which captures the names before the row and its AssetName children are
    deleted.

    The path is resolved here, while the row is still present, but nothing is
    removed until the transaction commits -- a rollback restores the row, and
    the content it points at must still be there.
    """
    if target.path is None:
        # Incomplete asset (no content on disk); nothing to remove.
        return
    genome_folder = _genome_folder(connection)
    content_dir = Path(genome_folder) / target.path
    if not content_dir.exists():
        # Warn, but let the deletion proceed: a catalog row for content that is
        # already gone is exactly the drift a delete is meant to resolve.
        logger.warning(
            f"Asset directory does not exist, but proceeding with deletion: {content_dir}"
        )
        return
    logger.info(f"Queued asset files for removal after commit: {content_dir}")
    enqueue_path_removal(target, content_dir, "tree")


def _alias_before_delete_handler(_, connection: Connection, target: Alias):
    """
    Queue the alias's on-disk trees for removal after the commit.

    Only fires for the SQL-backed AliasManager (server mode). Local mode stores
    aliases in the refget store, so StoreAliasManager queues the same paths
    itself — both routes share ``alias_owned_paths`` so they cannot diverge, and
    both remove files only after their catalog agrees the alias is gone.
    """
    genome_folder = _genome_folder(connection)
    for path in alias_owned_paths(genome_folder=Path(genome_folder), alias_name=target.name):
        enqueue_path_removal(target, path, "tree")


def _staged_asset_before_delete_handler(_, connection: Connection, target: StagedAsset):
    """
    Queue the staged asset's files for removal after the commit.
    Derives the directory from convention: genome_stage_folder / genome_digest / group_name.

    For mode="archive": removes the tarball (content-addressed by asset digest).
    For mode="file": removes the directory of per-file symlinks (originals in
    genome_folder untouched).
    """
    from refgenie.utils.staging import staged_archive_path

    genome_stage_folder = connection.execute(
        select(Configuration.genome_stage_folder).where(
            Configuration.id == select(func.max(Configuration.id)).scalar_subquery()
        )
    ).scalar_one_or_none()
    if not genome_stage_folder:
        logger.warning("No genome_stage_folder configured. Skipping file cleanup.")
        return

    asset_row = connection.execute(
        select(Asset.name, Asset.asset_group_id).where(Asset.digest == target.asset_digest)
    ).one_or_none()
    if not asset_row:
        logger.warning(f"Asset {target.asset_digest} not found. Skipping file cleanup.")
        return

    ag_row = connection.execute(
        select(AssetGroup.genome_digest, AssetGroup.name).where(
            AssetGroup.id == asset_row.asset_group_id
        )
    ).one_or_none()
    if not ag_row:
        return

    group_dir = Path(genome_stage_folder) / ag_row.genome_digest / ag_row.name

    if target.mode == "archive":
        # Tarball is content-addressed: {genome_digest}/{group}/{asset_digest}.tgz
        tarball_path = staged_archive_path(
            genome_stage_folder,
            ag_row.genome_digest,
            ag_row.name,
            target.asset_digest,
        )
        enqueue_path_removal(target, tarball_path, "file")
    elif target.mode == "file":
        # File mode stages group_name/asset_name/ as a directory of per-file
        # symlinks. (Legacy staged assets used a single directory symlink.)
        enqueue_path_removal(target, group_dir / asset_row.name, "tree")

    # Prune the parent group_dir, but only once the entries above are gone --
    # hence a callable queued behind them rather than a second path removal.
    enqueue_callable(target, lambda: _rmdir_if_empty(group_dir))


def _rmdir_if_empty(path: Path) -> None:
    if path.is_dir() and not any(path.iterdir()):
        path.rmdir()
        logger.info(f"Removed empty directory: {path}")


def _seek_key_before_insert_handler(_, connection: Connection, target: SeekKey):
    """
    Compute and set the seek key's on-disk size before insert; non-path seek
    keys get size=None.
    """
    from refgenie.db.tables import is_path_type

    if not is_path_type(target.type):
        # Non-path seek keys have no filesystem size
        target.size = None
        return

    genome_folder = _genome_folder(connection)
    seek_key_path = Path(genome_folder) / target.asset.path / target.value
    if target.type == SeekKeyType.file:
        size = seek_key_path.stat().st_size
    elif target.type == SeekKeyType.directory:
        size = sum(f.stat().st_size for f in seek_key_path.glob("**/*") if f.is_file())
    elif target.type == SeekKeyType.prefix:
        matched_files = (Path(genome_folder) / target.asset.path).glob(f"{target.value}*")
        size = sum(f.stat().st_size for f in matched_files if f.is_file())
    else:
        raise ValueError(f"Seek key size calculation not supported for {target.type=}")
    target.size = size
