"""
AssetPuller - handles pulling assets from remote refgenie servers.

This is an internal service class, not user-facing. ``rgc.transfer``
(``TransferManager``) owns it and uses it for ``transfer.pull`` and the bulk
pulls. It depends on ``rgc.asset``, ``rgc.servers`` and ``rgc.genome``; none of
them refer back to it.
"""

import shutil
import signal
import sys
import threading
from collections.abc import Callable
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING

from sqlalchemy.engine import Engine
from ubiquerg import untar  # type: ignore

from refgenie import progress
from refgenie.const import DEFAULT_PULL_SIZE_CUTOFF_GB
from refgenie.exceptions import (
    AssetExistsError,
    MissingAliasError,
    MissingAssetGroupError,
    NoArchiveError,
    PullFailedError,
    PullSkipped,
    ServerCannotServe,
)
from refgenie.logger import logger
from refgenie.managers.asset.colocation import recreate_colocation_symlinks
from refgenie.managers.base import ResourceManager
from refgenie.managers.genome import GenomeCreationResult
from refgenie.managers.sources.api_ids import (
    API_ID_ARCHIVE,
    API_ID_ASSET_RELATIONSHIPS,
)
from refgenie.models import GenomeAlias, GenomeDigest
from refgenie.plugins.events import NULL_EVENTS, EventSink
from refgenie.plugins.hooks import POST_PULL, PRE_PULL, HookEvent
from refgenie.utils.digest import BuildProvenance, checksum
from refgenie.utils.prompt import Confirmer, resolve_confirmer

if TYPE_CHECKING:
    from refgenie.db.tables import Asset
    from refgenie.managers.alias import AliasManager
    from refgenie.managers.asset.manager import AssetManager
    from refgenie.managers.genome import GenomeManager
    from refgenie.managers.sources.client import ServerClient
    from refgenie.managers.sources.servers import ServerManager


def _provenance_from_server(asset_names: list[dict], name: str) -> "BuildProvenance | None":
    """
    The build the server recorded for one of an asset's names.

    Without this a pulled asset's ``build_digest`` would exist only on the
    machine that built it, which defeats making builds addressable. Returns
    None when the server reported no build for that name.
    """
    entry = next((e for e in asset_names if e.get("name") == name), None)
    if entry is None or entry.get("build_digest") is None:
        return None
    timestamp = entry.get("build_timestamp")
    return BuildProvenance(
        build_digest=entry["build_digest"],
        build_level1=entry.get("build_level1"),
        build_digest_scheme=entry.get("build_digest_scheme"),
        build_timestamp=datetime.fromisoformat(timestamp) if timestamp else None,
        refgenie_version=entry.get("refgenie_version"),
        inputs=entry.get("inputs"),
        docker_image=entry.get("docker_image"),
        docker_image_digest=entry.get("docker_image_digest"),
        # Deliberately dropped: recipe ids are per-database surrogate keys and
        # the server's numbering says nothing about the client's.
        recipe_id=None,
    )


def handle_sigint_pull(filepath: Path):
    """Report an interrupted download and exit.

    Removing the partial file is not this handler's job: the download writes to
    a ``.part`` file beside ``filepath`` and unlinks it on any exception,
    including the ``SystemExit`` raised here. ``filepath`` itself never exists
    until the download has completed and been renamed into place.
    """

    def handle(sig, frame):
        logger.warning(f"\nThe download was interrupted: {filepath}")
        sys.exit(0)

    return handle


def should_pull_large_archive(
    archive_size: float | int,
    asset_registry_path: str,
    size_cutoff: float | int | None = None,
    force: bool | None = None,
    confirm: Confirmer | None = None,
) -> bool:
    """
    Check whether an archive over the size cutoff should be downloaded.

    Args:
        archive_size: The size of the archive, in bytes.
        asset_registry_path: The asset registry path.
        size_cutoff: The size cutoff, in GB.
        force: True to download regardless of size, False to skip, None to ask.
        confirm: Confirmation callback. Defaults to a refusal unless the CLI
            has enabled interactive prompts; see `refgenie.utils.prompt`.

    Returns:
        bool: Whether to pull the large archive.
    """
    size_cutoff_gb = size_cutoff if size_cutoff is not None else DEFAULT_PULL_SIZE_CUTOFF_GB
    size_cutoff_bytes = size_cutoff_gb * 1000**3
    logger.debug(f"'{asset_registry_path}' archive size: {archive_size}")
    if not force and archive_size > size_cutoff_bytes:
        if force is False:
            logger.info(f"Skipping pull of '{asset_registry_path}'; size: {archive_size}")
            return False
        # force is None, ask
        if not resolve_confirmer(confirm)(
            f"This archive exceeds the size cutoff "
            f"({archive_size / 1000**3:.1f}GB > {size_cutoff_gb:.1f}GB). Do you want to proceed?"
        ):
            logger.info(f"Skipping pull of '{asset_registry_path}'. Size too large: {archive_size}")
            return False
    return True


@dataclass
class PulledAssetMetadata:
    """Server metadata resolved for a single asset pull."""

    asset_metadata: dict
    asset_group_metadata: dict
    asset_digest: str
    declared_modes: list[str]
    asset_parents: list
    asset_name: str
    asset_class_name: str
    asset_names: list


@dataclass
class PullTransaction:
    """Context manager for atomic pull operations. Rolls back on failure."""

    alias_manager: "AliasManager"
    genome_manager: "GenomeManager"
    _created_genomes: list = field(default_factory=list)
    _created_aliases: list = field(default_factory=list)
    _created_dirs: list = field(default_factory=list)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is not None:
            self._rollback()
        return False

    def _rollback(self):
        """Roll back changes made during the pull operation."""
        logger.warning("Pull failed, rolling back changes...")
        # Remove directories in reverse order
        for dir_path in reversed(self._created_dirs):
            if dir_path.exists():
                shutil.rmtree(dir_path)
                logger.debug(f"Removed directory: {dir_path}")
        # Remove aliases we created
        for alias_name in reversed(self._created_aliases):
            try:
                self.alias_manager.remove(alias_name)
                logger.debug(f"Rolled back alias: {alias_name}")
            except Exception:
                pass  # Alias may already be gone
        # Remove genomes we created
        for genome_digest in reversed(self._created_genomes):
            try:
                self.genome_manager.remove(genome_digest)
                logger.debug(f"Rolled back genome: {genome_digest}")
            except Exception:
                pass  # Genome may already be gone

    def track_creation(self, result: GenomeCreationResult):
        """Accumulate genome/alias creations for potential rollback.

        pull() calls this once per alias when registering a multi-alias
        genome; rollback must undo all of them, not just the last one
        tracked.
        """
        if result.created_alias and result.alias_name:
            self._created_aliases.append(result.alias_name)
        if result.created_genome and result.genome_digest:
            self._created_genomes.append(result.genome_digest)

    def track_directory(self, path: Path):
        """Track a created directory for potential rollback."""
        self._created_dirs.append(path)


class AssetPuller(ResourceManager):
    """
    Handles pulling assets from remote refgenie servers.

    This is an internal service class used by TransferManager.pull(). It downloads
    assets from subscribed servers and records them through ``rgc.asset``.
    """

    def __init__(
        self,
        database_engine: Engine,
        asset: "AssetManager",
        servers: "ServerManager",
        alias_manager: "AliasManager",
        genome_manager: "GenomeManager",
        events: EventSink | None = None,
    ):
        """
        Initialize the AssetPuller.

        Args:
            database_engine: The database engine.
            asset: The AssetManager pulled assets are recorded in. Its
                folders, ``content``, ``group``, ``tree`` and ``links`` are
                what a pull writes through.
            servers: The ServerManager, for subscriptions and clients.
            alias_manager: The AliasManager for alias operations.
            genome_manager: The GenomeManager, which creates missing genomes.
            events: Where ``pre_pull`` / ``post_pull`` are emitted for plugins.
        """
        super().__init__(database_engine)
        self._events = events or NULL_EVENTS
        self._asset = asset
        self._servers = servers
        self._alias_manager = alias_manager
        self._genome_manager = genome_manager

    def _alias_tree_is_missing(
        self,
        alias_name: GenomeAlias,
        asset_group_name: str,
        asset_name: str,
    ) -> bool:
        """
        Whether the name-addressed view of an asset is absent for this alias.

        The alias tree is derived from the catalog, so its absence is never
        information -- it only ever means the rendering step did not happen. That
        makes it safe to re-render, and makes refusing a pull because the row
        exists the wrong answer when the tree does not.

        Args:
            alias_name: The genome alias the pull was addressed to.
            asset_group_name: The asset group name.
            asset_name: The asset name.

        Returns:
            bool: True when the alias directory for this asset does not exist.
        """
        from refgenie.utils.symlinks import get_symlink_paths

        paths = get_symlink_paths(
            alias_folder=self._asset.alias_folder,
            aliases=[alias_name],
            asset_group_name=asset_group_name,
            asset_name=asset_name,
        )
        return not all(path.exists() for path in paths.values())

    # === Bulk Pull Methods ===

    def pull_multiple(
        self,
        asset_list: list[dict],
        force: bool | None = None,
        force_large: bool | None = None,
        size_cutoff: int | float | None = None,
        sigint_handler: Callable | None = None,
        confirm: Confirmer | None = None,
    ) -> list:
        """
        Pull multiple assets. Calls self.pull() for each.

        Continues on individual failures, reports summary at end.

        Args:
            asset_list: Asset dicts from ServerManager.list_assets_for_genome().
                        Each dict should have: genome_digest, asset_group_name, asset_name.
            force: How to handle case in which asset path already exists.
            force_large: How to handle case of large archive.
            size_cutoff: Maximum archive file size to download without prompt.
            sigint_handler: Signal handler for interrupts during download.

        Returns:
            List of successfully pulled Assets.
        """
        from refgenie.db.tables import Asset

        successful: list[Asset] = []
        failed: list[dict] = []

        for asset_info in asset_list:
            genome_digest = asset_info.get("genome_digest")
            asset_group_name = asset_info.get("asset_group_name")
            asset_name = asset_info.get("asset_name")
            server_url = asset_info.get("server_url")

            if not all([genome_digest, asset_group_name]):
                logger.warning(f"Skipping invalid asset entry: {asset_info}")
                failed.append(asset_info)
                continue

            logger.info(f"Pulling {genome_digest}/{asset_group_name}:{asset_name or 'default'}")
            try:
                result = self.pull(
                    asset_group_name=asset_group_name,
                    genome=GenomeDigest(genome_digest),
                    asset_name=asset_name,
                    force=force,
                    force_large=force_large,
                    force_server_urls=[server_url] if server_url else None,
                    size_cutoff=size_cutoff,
                    sigint_handler=sigint_handler,
                    confirm=confirm,
                )
                if result is not None:
                    successful.append(result)
                else:
                    failed.append(asset_info)
            except Exception as e:
                logger.error(f"Failed to pull {genome_digest}/{asset_group_name}:{asset_name}: {e}")
                failed.append(asset_info)

        # Report summary
        if failed:
            logger.warning(f"Pull summary: {len(successful)} succeeded, {len(failed)} failed")
        else:
            logger.info(f"Pull summary: {len(successful)} assets pulled successfully")

        return successful

    # === Pull Methods ===

    def pull(
        self,
        asset_group_name: str,
        genome: GenomeAlias | GenomeDigest,
        asset_name: str | None = None,
        force: bool | None = None,
        force_large: bool | None = None,
        force_server_urls: list[str] | None = None,
        size_cutoff: int | float | None = None,
        sigint_handler: Callable | None = None,
        confirm: Confirmer | None = None,
    ):
        """
        Download and register one asset, firing ``pre_pull`` and ``post_pull``.

        ``post_pull`` fires on every outcome: ``succeeded`` is True when an
        asset came back or was already present, and False when the pull
        returned nothing or raised. See :meth:`_pull` for the pull itself.

        Raises:
            TypeError: If ``genome`` is a plain ``str`` rather than one of the
                two genome identifier types.
        """
        if not isinstance(genome, (GenomeAlias, GenomeDigest)):
            raise TypeError(
                f"pull needs a GenomeAlias or a GenomeDigest, not {type(genome).__name__} "
                f"'{genome}'. Wrap it to say which it is."
            )
        self._events.emit(
            HookEvent(
                hook=PRE_PULL, genome=str(genome), asset_group=asset_group_name, asset=asset_name
            )
        )
        result = None
        succeeded = False
        try:
            result = self._pull(
                asset_group_name=asset_group_name,
                genome=genome,
                asset_name=asset_name,
                force=force,
                force_large=force_large,
                force_server_urls=force_server_urls,
                size_cutoff=size_cutoff,
                sigint_handler=sigint_handler,
                confirm=confirm,
            )
            succeeded = result is not None
            return result
        except AssetExistsError:
            succeeded = True  # the asset is here; nothing needed pulling
            raise
        finally:
            self._events.emit(
                HookEvent(
                    hook=POST_PULL,
                    genome=str(genome),
                    asset_group=asset_group_name,
                    asset=getattr(result, "name", None) or asset_name,
                    succeeded=succeeded,
                )
            )

    def _pull(
        self,
        asset_group_name: str,
        genome: GenomeAlias | GenomeDigest,
        asset_name: str | None = None,
        force: bool | None = None,
        force_large: bool | None = None,
        force_server_urls: list[str] | None = None,
        size_cutoff: int | float | None = None,
        sigint_handler: Callable | None = None,
        confirm: Confirmer | None = None,
    ):
        """
        Download and register one asset from the first subscribed server that
        can serve it.

        Uses injected managers for all cross-manager operations.
        Uses PullTransaction for atomic rollback on failure.

        Args:
            asset_group_name: Name of a group of assets to fetch.
            genome: The genome, as a ``GenomeAlias`` or a ``GenomeDigest``. It
                need not exist locally; see :meth:`_resolve_pull_genome`.
            asset_name: Name of particular asset to fetch.
            force: How to handle case in which asset path already exists.
            force_large: How to handle archives larger than size_cutoff (default 10GB).
            force_server_urls: Force specific server URLs to use.
            size_cutoff: Maximum archive file size to download without prompt.
            sigint_handler: Signal handler for interrupts during download.

        Returns:
            Asset or None: The added asset, or None if pull failed.
        """
        added_asset = None
        last_error = None

        # The caller (TransferManager.pull) handles the subscribe prompt.
        server_urls = force_server_urls or self._servers.subscriptions()
        if not server_urls:
            raise PullFailedError("No server subscriptions found")

        for server_url in server_urls:
            # Each attempt starts from the caller's request: the attempt's
            # server-specific resolutions (canonical asset name, digest, alias)
            # are locals of `_pull_from_server`, so a failed attempt cannot leak
            # them into the next one.
            try:
                added_asset = self._pull_from_server(
                    server_url=server_url,
                    asset_group_name=asset_group_name,
                    genome=genome,
                    asset_name=asset_name,
                    force=force,
                    force_large=force_large,
                    size_cutoff=size_cutoff,
                    sigint_handler=sigint_handler,
                    confirm=confirm,
                )
            except ServerCannotServe as exc:
                # The transaction has already rolled back whatever this attempt
                # created; move on to the next server.
                logger.warning(f"Server {server_url} cannot serve this asset: {exc}")
                if exc.final_error is not None:
                    last_error = exc.final_error
                continue
            except PullSkipped as exc:
                # Deliberate skip -- do not try other servers.
                logger.info(str(exc))
                return None

            # This server served the asset. Stop here: re-entering the loop
            # would re-check a now-registered asset and raise AssetExistsError
            # for a pull that just succeeded.
            break

        if added_asset is None and last_error:
            raise last_error
        return added_asset

    def _pull_from_server(
        self,
        server_url: str,
        asset_group_name: str,
        genome: GenomeAlias | GenomeDigest,
        asset_name: str | None,
        force: bool | None,
        force_large: bool | None,
        size_cutoff: int | float | None,
        sigint_handler: Callable | None,
        confirm: Confirmer | None,
    ) -> Asset | None:
        """
        One pull attempt against one server, inside its own transaction.

        Raises:
            ServerCannotServe: This server cannot serve the asset; the
                transaction has rolled back and the caller may try the next one.
            PullSkipped: The user or ``force=False`` declined to replace
                existing content.
        """
        # Use PullTransaction for atomic rollback on failure
        with PullTransaction(self._alias_manager, self._genome_manager) as txn:
            client = self._servers.client(server_url)
            alias_name, genome_digest = self._resolve_pull_genome(genome, client, server_url, txn)
            try:
                asset_name = str(
                    asset_name
                    or self._asset.group.get_default(asset_group_name, genome_digest=genome_digest)
                )
            except MissingAssetGroupError:
                pass
            else:
                if not force and self._asset.exists(
                    genome_digest=genome_digest,
                    asset_group_name=asset_group_name,
                    asset_name=asset_name,
                ):
                    # The row is committed but the alias tree is not
                    # rendered: a pull killed between the two. The tree
                    # is derived from the catalog, so render it and call
                    # the pull done.
                    if self._alias_tree_is_missing(
                        alias_name=alias_name,
                        asset_group_name=asset_group_name,
                        asset_name=asset_name,
                    ):
                        logger.warning(
                            f"Asset '{alias_name}/{asset_group_name}:{asset_name}' is in "
                            f"the catalog but its alias tree is missing. Rendering it."
                        )
                        self._asset.tree.render(genome_digest, asset_group_name, asset_name)
                        return self._asset.get(
                            genome_digest=genome_digest,
                            asset_group_name=asset_group_name,
                            asset_name=asset_name,
                        )
                    raise AssetExistsError(
                        f"Asset '{alias_name}/{asset_group_name}:{asset_name}' already exists"
                    )

            if asset_name is None:
                logger.debug(
                    f"No local default asset for '{asset_group_name}'; "
                    f"will accept any single asset from the server"
                )

            meta = self._fetch_asset_metadata(
                client=client,
                server_url=server_url,
                genome_digest=genome_digest,
                asset_group_name=asset_group_name,
                asset_name=asset_name,
                alias_name=alias_name,
            )
            # The server-metadata section resolves the canonical asset name;
            # adopt it for everything downstream.
            asset_name = meta.asset_name
            asset_digest = meta.asset_digest
            asset_parents = meta.asset_parents

            bundle_name = f"{alias_name}/{asset_group_name}:{asset_name}"

            # Content is addressed by digest, both here and on the server, so
            # the download directory is the digest directory. The archive's
            # top-level entry is that same digest, so extraction lands here.
            asset_dir = self._asset.data_folder / genome_digest / asset_group_name / asset_digest
            self._clear_existing_content(asset_dir, asset_digest, bundle_name, force, confirm)

            # Remove existing database record when force-overwriting.
            # Done after directory cleanup (above) to avoid the ORM event
            # handler deleting newly downloaded files.
            if force and self._asset.exists(
                genome_digest=genome_digest,
                asset_group_name=asset_group_name,
                asset_name=asset_name,
            ):
                self._asset.remove(
                    genome_digest=genome_digest,
                    asset_group_name=asset_group_name,
                    asset_name=asset_name,
                )

            if asset_parents:
                # Check if parent digests exist locally
                if not self._asset.links.parent_digests_exist(
                    asset_parents,
                    lambda digest: self._asset.get_by_digest(digest) is not None,
                ):
                    raise PullFailedError("Parent assets not available locally")

            self._download_staged(
                client=client,
                server_url=server_url,
                asset_digest=asset_digest,
                declared_modes=meta.declared_modes,
                genome_digest=genome_digest,
                asset_group_name=asset_group_name,
                asset_name=asset_name,
                asset_dir=asset_dir,
                bundle_name=bundle_name,
                txn=txn,
                force_large=force_large,
                size_cutoff=size_cutoff,
                sigint_handler=sigint_handler,
                confirm=confirm,
            )

            # Post-download registration; must stay inside the transaction
            # so rollback covers it.
            return self._finalize_pulled_asset(
                client=client,
                alias_name=alias_name,
                genome_digest=genome_digest,
                asset_group_name=asset_group_name,
                asset_name=asset_name,
                asset_dir=asset_dir,
                meta=meta,
            )

    def _clear_existing_content(
        self,
        asset_dir: Path,
        asset_digest: str,
        bundle_name: str,
        force: bool | None,
        confirm: Confirmer | None,
    ) -> None:
        """
        Make way for a download into ``asset_dir``, or raise ``PullSkipped``.

        Orphaned content (no complete catalog row) is always reclaimed; real
        content is replaced only when ``force`` or the confirmer allows it.
        """
        if asset_dir.exists():
            logger.info(f"Asset directory already exists: {asset_dir}")
            existing_row = self._asset.get_by_digest(asset_digest)
            if existing_row is None or existing_row.path is None:
                # Content on disk with no complete catalog row is not
                # an asset, it is debris from a pull that was killed
                # between extraction and the commit. Reclaim it.
                logger.warning(
                    f"Reclaiming an orphaned asset directory with no catalog entry: {asset_dir}"
                )
                shutil.rmtree(asset_dir)
            elif force is False:
                raise PullSkipped(f"Skipping pull of {bundle_name}")
            elif force is None:
                # Go through the caller's confirmer, never a bare
                # rich prompt: a pull driven from a server request
                # or a worker thread has no terminal to answer
                # with, and `Confirm.ask` there blocks forever.
                # `resolve_confirmer` returns `deny` unless the CLI
                # enabled interactive prompts.
                if not resolve_confirmer(confirm)(f"Replace existing {asset_dir}?"):
                    raise PullSkipped(f"Skipping pull of '{bundle_name}'")
                else:
                    logger.debug(f"Overwriting: {asset_dir}")
                    shutil.rmtree(asset_dir)
            else:
                logger.debug(f"Overwriting: {asset_dir}")
                shutil.rmtree(asset_dir)

    def _download_staged(
        self,
        client: "ServerClient",
        server_url: str,
        asset_digest: str,
        declared_modes: list[str],
        genome_digest: GenomeDigest,
        asset_group_name: str,
        asset_name: str,
        asset_dir: Path,
        bundle_name: str,
        txn: PullTransaction,
        force_large: bool | None,
        size_cutoff: int | float | None,
        sigint_handler: Callable | None,
        confirm: Confirmer | None,
    ) -> None:
        """
        Download whatever the server has staged for ``asset_digest`` into ``asset_dir``.

        Raises:
            ServerCannotServe: Nothing downloadable is staged on this server.
            PullFailedError: The asset is metadata-only, or staging is ambiguous.
        """
        # Inspect what is actually staged for this asset. A single,
        # unfiltered query lets us branch on real availability instead of
        # trusting a mode filter that may not match what was staged.
        try:
            staged_items = client.get_staged_assets(params={"asset_digest": asset_digest})
        except Exception as e:
            raise ServerCannotServe(f"Failed to query staged assets from {server_url}: {e}") from e

        archive_records = [rec for rec in staged_items if rec["mode"] == "archive"]
        file_records = [rec for rec in staged_items if rec["mode"] == "file"]

        # Branch on what is downloadable.
        if archive_records and "archive" in declared_modes:
            if len(archive_records) != 1:
                raise PullFailedError(
                    f"Expected one archive, got {len(archive_records)} for "
                    f"asset digest {asset_digest} on server {server_url}"
                )
            self._pull_archive_mode(
                client=client,
                staged_asset_metadata=archive_records[0],
                asset_digest=asset_digest,
                genome_digest=genome_digest,
                asset_group_name=asset_group_name,
                asset_name=asset_name,
                asset_dir=asset_dir,
                bundle_name=bundle_name,
                txn=txn,
                force_large=force_large,
                size_cutoff=size_cutoff,
                sigint_handler=sigint_handler,
                confirm=confirm,
            )
        elif file_records:
            self._pull_file_mode(client, asset_digest, asset_dir, bundle_name, txn)
        elif declared_modes == ["none"]:
            raise PullFailedError(
                f"Asset '{bundle_name}' is metadata-only (serving_modes=['none']). "
                f"Use 'refgenie build' to create it locally."
            )
        else:
            message = (
                f"No archive found for asset digest {asset_digest} "
                f"on server {server_url}. "
                f"The asset may not be available for download."
            )
            raise ServerCannotServe(message, final_error=NoArchiveError(message))

    def _resolve_pull_genome(
        self,
        genome: GenomeAlias | GenomeDigest,
        client: "ServerClient",
        server_url: str,
        txn: PullTransaction,
    ) -> tuple[GenomeAlias, GenomeDigest]:
        """
        The alias and the digest a pull of ``genome`` is addressed to.

        A pull needs an alias as well as a digest, because the alias tree it
        renders is keyed by one. The type of ``genome`` says which one the
        caller has:

        - a ``GenomeAlias`` is used as given. If it is not known here, the
          server resolves it and it is registered locally, which is how a
          first pull creates the genome.
        - a ``GenomeDigest`` is registered under the aliases this server
          reports for it, and the pull goes to the first one (else to a local
          alias).

        Everything registered here is tracked on ``txn``, so a failed pull
        rolls it back.

        Raises:
            ServerCannotServe: If this server cannot identify the genome.
            PullFailedError: If the server's alias listing fails.
        """
        if isinstance(genome, GenomeAlias):
            return self._resolve_pull_alias(genome, server_url, txn)
        return self._resolve_pull_digest(genome, client, server_url, txn)

    def _resolve_pull_alias(
        self,
        genome_alias: GenomeAlias,
        server_url: str,
        txn: PullTransaction,
    ) -> tuple[GenomeAlias, GenomeDigest]:
        """The pull target for an alias: a local one, else one the server knows."""
        try:
            return genome_alias, self._alias_manager.resolve(genome_alias)
        except MissingAliasError:
            pass
        logger.warning(
            f"No local genome alias '{genome_alias}'. "
            f"Setting genome identity with server: {server_url}"
        )
        result = self._genome_manager.ensure_from_remote(
            alias_name=genome_alias, server_urls=[server_url]
        )
        if not result.success:
            raise ServerCannotServe(
                f"Could not resolve genome alias '{genome_alias}' via server {server_url}."
            )
        txn.track_creation(result)
        return genome_alias, self._alias_manager.resolve(genome_alias)

    def _resolve_pull_digest(
        self,
        genome_digest: GenomeDigest,
        client: "ServerClient",
        server_url: str,
        txn: PullTransaction,
    ) -> tuple[GenomeAlias, GenomeDigest]:
        """The pull target for a digest: the server's aliases for it, else a local one."""
        known_digest = self._genome_manager.exists(genome_digest)
        try:
            # Use get_all_aliases method (required by ServerClient protocol)
            aliases_metadata = client.get_all_aliases(params={"genome_digest": genome_digest})
        except Exception as e:
            if not known_digest:
                raise ServerCannotServe(
                    f"Could not resolve genome digest '{genome_digest}' via server {server_url}."
                ) from e
            raise PullFailedError(f"Failed to fetch aliases metadata from {server_url}: {e}") from e
        for alias_metadata in aliases_metadata:
            result = self._genome_manager.ensure_from_remote(
                alias_name=GenomeAlias(str(alias_metadata["name"])),
                genome_digest=GenomeDigest(alias_metadata["genome_digest"]),
                genome_description=alias_metadata.get("genome_description"),
            )
            txn.track_creation(result)
        if aliases_metadata:
            # Address the pull to the first alias, not whichever entry happened
            # to be listed last.
            return GenomeAlias(str(aliases_metadata[0]["name"])), genome_digest
        local_aliases = self._alias_manager.get_for_genome(genome_digest) if known_digest else []
        if local_aliases:
            return GenomeAlias(local_aliases[0]), genome_digest
        raise ServerCannotServe(
            f"Could not resolve genome digest '{genome_digest}' via server {server_url}."
        )

    def _fetch_asset_metadata(
        self,
        client: "ServerClient",
        server_url: str,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        asset_name: str | None,
        alias_name: GenomeAlias | None,
    ) -> PulledAssetMetadata:
        """
        Query the server for asset group, asset, and relationship metadata.

        Resolves the canonical asset name from server metadata (returned in the
        bundle rather than mutating the caller's variable).
        """
        logger.info(
            f"Querying server {server_url} for {genome_digest}/{asset_group_name}"
            f"{':' + asset_name if asset_name else ''}"
        )

        # Query /asset_groups to get asset group ID
        try:
            asset_groups_items = client.get_asset_groups(
                params={
                    "genome_digest": genome_digest,
                    "asset_group_name": asset_group_name,
                }
            )
        except Exception as e:
            raise ServerCannotServe(f"Failed to query asset groups from {server_url}: {e}") from e

        if len(asset_groups_items) != 1:
            raise PullFailedError(
                f"Expected one asset group, got {len(asset_groups_items)} for "
                f"{genome_digest}/{asset_group_name} on server {server_url}"
            )
        asset_group_metadata = asset_groups_items[0]
        asset_group_id = asset_group_metadata["id"]

        # Query /assets to get asset metadata using asset_group_id
        try:
            assets_items = client.get_assets(
                params={
                    "asset_group_id": asset_group_id,
                    "name": asset_name,
                }
            )
        except Exception as e:
            raise ServerCannotServe(f"Failed to query assets from {server_url}: {e}") from e

        if len(assets_items) != 1:
            if asset_name is None and len(assets_items) > 1:
                names = [a.get("name", "?") for a in assets_items]
                raise PullFailedError(
                    f"Multiple assets found for {genome_digest}/{asset_group_name}: {names}. "
                    f"Specify one with: refgenie pull {alias_name}/{asset_group_name}:<asset_name>"
                )
            else:
                raise PullFailedError(
                    f"Expected one asset, got {len(assets_items)} for "
                    f"asset_group_id {asset_group_id}, asset_name {asset_name} "
                    f"on server {server_url}"
                )

        asset_metadata = assets_items[0]
        asset_digest = asset_metadata["digest"]
        # The server is the single source of truth for the resolved serving
        # modes; the client does not re-derive them from serving_modes_override.
        declared_modes = asset_metadata["serving_modes"]

        # Get relationship metadata for parent assets
        relationship_metadata = client.get(
            operation_id=API_ID_ASSET_RELATIONSHIPS,
            url_format_params={"asset_digest": asset_digest},
            params={"expand": True},
        )
        asset_parents = relationship_metadata["parents"]
        asset_name = str(asset_metadata["name"])
        asset_class_name = asset_group_metadata["name"]
        # Every name the server knows this content by, each with the build
        # provenance behind it (AssetNameResponse).
        asset_names = asset_metadata["names"]

        return PulledAssetMetadata(
            asset_metadata=asset_metadata,
            asset_group_metadata=asset_group_metadata,
            asset_digest=asset_digest,
            declared_modes=declared_modes,
            asset_parents=asset_parents,
            asset_name=asset_name,
            asset_class_name=asset_class_name,
            asset_names=asset_names,
        )

    def _finalize_pulled_asset(
        self,
        client: "ServerClient",
        alias_name: GenomeAlias,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        asset_name: str,
        asset_dir: Path,
        meta: PulledAssetMetadata,
    ):
        """
        Register a downloaded asset: colocation symlinks, fasta genome
        initialization, catalog write (renders the alias tree), extra-name
        adoption, default-asset set, and parent relationships.

        Must be called inside the PullTransaction context so rollback covers it.

        Returns:
            Asset or None: The added asset.
        """
        progress.emit("stage", message="Registering asset", phase="register")
        asset_metadata = meta.asset_metadata
        asset_digest = meta.asset_digest
        asset_class_name = meta.asset_class_name
        asset_names = meta.asset_names
        asset_parents = meta.asset_parents

        # Post-download: recreate colocation symlinks
        colocate_metadata = asset_metadata.get("colocate")
        if colocate_metadata:
            recreate_colocation_symlinks(
                output_folder=asset_dir,
                genome_folder=self._asset.genome_folder,
                colocate_metadata=colocate_metadata,
                parent_assets=self._resolve_colocation_parents(colocate_metadata, genome_digest),
            )

        # Post-download: fasta genome initialization
        if asset_class_name == "fasta":
            local_digest, _created = self._genome_manager.initialize_genome(
                fasta_file_path=asset_dir / f"{genome_digest}.fa",
                # Used when the genome is (re)created here, e.g. after a forced
                # re-pull removed the old fasta and, with it, the genome.
                description=self._servers.genome_description(genome_digest, [client.server_url]),
                alias_names=[alias_name],
                use_existing=True,
            )
            if local_digest != genome_digest:
                raise ValueError(
                    f"Digest mismatch for '{alias_name}': "
                    f"server={genome_digest}, local={local_digest}. "
                    f"Your refgenie version computes a different digest than the server. "
                    f"To fix:\n"
                    f"  pip install -U refgenie\n"
                    f"  refgenie remove {alias_name}/fasta -f\n"
                    f"  refgenie pull {alias_name}/fasta"
                )

        # Add the downloaded asset under the server's digest rather than
        # recomputing it. The building recipe's `inherent` set
        # decides which files the digest covers, and that declaration
        # does not travel with the archive -- so a client recomputing it
        # would disagree with the server for any recipe that declares
        # one. `asset_digest` is already the authority for the content
        # directory above; this makes it the authority for the identity
        # too. Content integrity is checked separately, against the
        # archive digest during download. content.add renders the alias tree
        # itself, after its commit.
        progress.emit("stage", message="Registering and linking", phase="symlink")
        added_asset = self._asset.content.add(
            asset_class_name=asset_class_name,
            path=asset_dir.relative_to(self._asset.genome_folder),
            genome_digest=genome_digest,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
            colocate=colocate_metadata,
            digest=asset_digest,
            build_provenance=_provenance_from_server(asset_names, asset_name),
        )

        # Adopt every other name the server knows this content by, so
        # `:0.7.19` resolves locally even when the server's canonical
        # name differs. A conflicting name is skipped with a warning, not
        # a failed pull.
        for entry in asset_names:
            extra_name = str(entry["name"])
            if extra_name == asset_name:
                continue
            self._asset.content.adopt_name(
                genome_digest=genome_digest,
                asset_group_name=asset_group_name,
                asset_name=extra_name,
                asset_digest=asset_digest,
                build_provenance=_provenance_from_server(asset_names, extra_name),
            )

        # Set default asset if none exists
        try:
            default_asset = self._asset.group.get_default(
                asset_group_name, genome_digest=genome_digest
            )
        except MissingAssetGroupError:
            default_asset = None
        if default_asset is None:
            self._asset.group.set_default(asset_group_name, asset_name, genome_digest=genome_digest)

        # Set relationships between the asset and its parents
        if asset_parents:
            self._asset.links.set_parents(
                genome_digest=genome_digest,
                asset_group_name=asset_group_name,
                asset_name=asset_name,
                parent_asset_digests=[p["digest"] for p in asset_parents],
            )

        return added_asset

    def _pull_archive_mode(
        self,
        client: "ServerClient",
        staged_asset_metadata: dict,
        asset_digest: str,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        asset_name: str,
        asset_dir: Path,
        bundle_name: str,
        txn: PullTransaction,
        force_large: bool | None,
        size_cutoff: int | float | None,
        sigint_handler: Callable | None,
        confirm: Confirmer | None = None,
    ) -> None:
        """Pull an asset by downloading its archive tarball.

        Args:
            client: The server client.
            staged_asset_metadata: Metadata dict for the StagedAsset with mode='archive'.
            asset_digest: The digest of the asset.
            genome_digest: The digest of the genome.
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            asset_dir: Local directory where asset files will be extracted.
            bundle_name: Human-readable asset identifier for logging.
            txn: The PullTransaction for rollback tracking.
            force_large: How to handle large archives.
            size_cutoff: Maximum archive file size to download without prompt.
            sigint_handler: Signal handler for interrupts during download.
        """
        if not should_pull_large_archive(
            archive_size=staged_asset_metadata.get("tarball_size"),
            asset_registry_path=bundle_name,
            force=force_large,
            size_cutoff=size_cutoff,
            confirm=confirm,
        ):
            raise RuntimeError(
                f"Archive download for '{bundle_name}' was skipped by user or size policy."
            )

        tardir = asset_dir.parent
        # Download the tarball into the group directory (a sibling of the
        # extracted content dir), not into asset_dir itself: the archive's
        # top-level entry IS asset_dir's basename (the content digest), so
        # untarring into tardir recreates asset_dir.
        tarpath = tardir / f"{asset_group_name}__{asset_name}.tgz"

        tardir_existed = tardir.exists()
        if not tardir_existed:
            logger.debug(f"Creating directory: {tardir}")
            tardir.mkdir(parents=True, exist_ok=True)
            txn.track_directory(tardir)

        # Download the archive. signal.signal() only works on the main thread
        # of the main interpreter (see the matching guard/comment in
        # builder.py); off the main thread, skip installation and leave
        # whatever handler is already there. Restore the previous handler on
        # the way out -- success or exception -- so the process's SIGINT does
        # not stay permanently bound to this download's tarball path.
        install_sigint = threading.current_thread() is threading.main_thread()
        previous_sigint_handler = None
        if install_sigint:
            sigint_handler = sigint_handler or handle_sigint_pull
            previous_sigint_handler = signal.signal(signal.SIGINT, sigint_handler(tarpath))
        try:
            client.download_with_progress(
                operation_id=API_ID_ARCHIVE,
                output_path=tarpath,
                url_format_params={"asset_digest": asset_digest},
                name=bundle_name,
            )
        finally:
            if install_sigint:
                signal.signal(signal.SIGINT, previous_sigint_handler)
        # Verify checksum. `tarball_digest` is the digest of the archive's BYTES,
        # not the asset identity digest (`digest`) -- the two are different values
        # and must never be substituted for one another.
        # A full sha256 pass over a multi-gigabyte archive, emitting nothing
        # while it runs. Announce it, or the UI looks hung.
        progress.emit("stage", message=f"Verifying {bundle_name}", phase="verify")
        server_digest = staged_asset_metadata.get("tarball_digest")
        if server_digest is None:
            logger.warning(
                f"Archive for '{bundle_name}' could not be verified: the server's "
                f"staged-asset record carries no 'tarball_digest'."
            )
        elif (local_digest := checksum(tarpath)) != server_digest:
            raise ValueError(
                f"Downloaded archive ({tarpath}) checksum mismatch: "
                f"({local_digest}, {server_digest})"
            )
        # Extract tarball
        if tarpath.suffix == ".tgz":
            progress.emit("stage", message=f"Unpacking {bundle_name}", phase="unpack")
            logger.info(f"Extracting asset tarball: {tarpath}")
            untar(tarpath.as_posix(), tardir.as_posix(), filter="fully_trusted")
            tarpath.unlink()

    def _pull_file_mode(
        self,
        client: "ServerClient",
        asset_digest: str,
        asset_dir: Path,
        bundle_name: str,
        txn: PullTransaction,
    ) -> None:
        """Pull an asset by downloading individual files from the file-level endpoint.

        Used when the asset's serving modes include 'file' but not 'archive'.

        Args:
            client: The server client.
            asset_digest: The digest of the asset.
            asset_dir: Local directory to download files into.
            bundle_name: Human-readable asset identifier for logging.
            txn: The PullTransaction for rollback tracking.
        """
        file_list = client.get_asset_file_list(asset_digest)
        if not file_list:
            msg = f"No files found for asset {asset_digest}"
            logger.error(msg)
            raise RuntimeError(msg)

        logger.info(f"Downloading {len(file_list)} file(s) for '{bundle_name}'")

        # Create asset directory and track for rollback
        if not asset_dir.exists():
            asset_dir.mkdir(parents=True, exist_ok=True)
            txn.track_directory(asset_dir)

        for file_path in file_list:
            output_path = asset_dir / file_path
            # Create parent directories within asset_dir as needed
            output_path.parent.mkdir(parents=True, exist_ok=True)
            client.download_file(asset_digest, file_path, output_path)

        logger.info(f"Downloaded {len(file_list)} file(s) for '{bundle_name}' to {asset_dir}")

    def _resolve_colocation_parents(
        self,
        colocate_metadata: list[dict[str, str]],
        genome_digest: GenomeDigest,
    ) -> dict[str, "Asset | None"]:
        """Resolve each colocation entry's parent asset group to its default asset.

        Args:
            colocate_metadata: Entries with a ``parent_asset_group``.
            genome_digest: The genome the parents belong to.

        Returns:
            Parent asset group name -> the resolved asset, or None if the group
            has no asset here (colocation for it is then skipped, with a warning).
        """
        parents: dict[str, Asset | None] = {}
        for entry in colocate_metadata:
            parent_group = entry.get("parent_asset_group")
            if not parent_group or parent_group in parents:
                continue
            try:
                parents[parent_group] = self._asset.get(
                    genome_digest=genome_digest,
                    asset_group_name=parent_group,
                    asset_name=self._asset.group.get_default(
                        asset_group_name=parent_group,
                        genome_digest=genome_digest,
                    )
                    or "default",
                )
            except Exception as e:
                logger.warning(f"Could not find parent asset '{parent_group}': {e}")
                parents[parent_group] = None
        return parents
