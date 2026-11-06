"""
AssetPuller - handles pulling assets from remote refgenie servers.

This is an internal service class, not user-facing. ``rgc.transfer``
(``TransferManager``) owns it and uses it for ``transfer.pull`` and the bulk
pulls. It depends on ``rgc.asset``, ``rgc.servers`` and ``rgc.genome``; none of
them refer back to it.
"""

import shutil
from dataclasses import dataclass, field, replace
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING

from refgenie.db.engine import Engine
from ubiquerg import untar  # type: ignore

from refgenie import progress
from refgenie.const import DEFAULT_PULL_SIZE_CUTOFF_GB
from refgenie.db.tables import SeekKeyType
from refgenie.exceptions import (
    AssetExistsError,
    FastaDigestMismatchError,
    MissingAliasError,
    MissingAssetClassError,
    MissingAssetError,
    MissingAssetGroupError,
    NoArchiveError,
    NoSubscriptionsError,
    ParentPullError,
    PullFailedError,
    PullSkipped,
    ServerCannotServe,
)
from refgenie.logger import logger
from refgenie.managers.asset.colocation import recreate_colocation_symlinks
from refgenie.managers.base import ResourceManager
from refgenie.managers.genome import GenomeCreationResult
from refgenie.managers.sources.api_ids import API_ID_ARCHIVE
from refgenie.models import GenomeAlias, GenomeDigest
from refgenie.plugins.events import NULL_EVENTS, EventSink
from refgenie.plugins.hooks import POST_PULL, PRE_PULL, HookEvent
from refgenie.utils.digest import BuildProvenance, checksum
from refgenie.utils.prompt import Confirmer, resolve_confirmer
from refgenie.utils.signals import sigint_guard

if TYPE_CHECKING:
    from refgenie.db.tables import Asset
    from refgenie.managers.alias import AliasManager
    from refgenie.managers.asset.manager import AssetManager
    from refgenie.managers.asset_class import AssetClassManager
    from refgenie.managers.genome import GenomeManager
    from refgenie.managers.sources.client import RefgenieserverClient
    from refgenie.managers.sources.servers import ServerManager


def _ask(server_url: str, what: str, fn, *args, **kwargs):
    """Call the server; any failure becomes ServerCannotServe so the next server is tried."""
    try:
        return fn(*args, **kwargs)
    except Exception as e:
        raise ServerCannotServe(f"Failed to query {what} from {server_url}: {e}") from e


def _provenance_from_server(asset_names: list[dict], name: str) -> "BuildProvenance | None":
    """The build the server recorded for one of an asset's names, or None."""
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
        recipe_id=None,  # server's recipe ids mean nothing in this database
    )


def _seek_keys_from_server(seek_keys: list[dict] | None) -> dict[str, tuple[str, SeekKeyType]]:
    """
    The server's seek key values for an asset, in ``custom_seek_keys`` form.

    Non-path values come from the building recipe, which the client never sees,
    so they must be taken from the server rather than re-derived.
    """
    return {sk["name"]: (sk["value"], SeekKeyType(sk["type"])) for sk in seek_keys or []}


def _unresolved_alias_reason(
    alias: GenomeAlias, client: "RefgenieserverClient", server_url: str
) -> str:
    """
    Why a server did not resolve ``alias``, in words the user can act on.

    Only a 404 is the server's answer. The alias lookup returns 404 both for an
    unknown alias and for one whose sequence collection the server cannot
    return; the alias listing tells them apart. Any other error is reported as is.
    """
    try:
        client.alias.resolve(alias)
    except Exception as e:  # noqa: BLE001 - reported to the user
        # By status, not class: an injected test client raises httpx2's error.
        if getattr(getattr(e, "response", None), "status_code", None) != 404:
            return f"Could not resolve genome alias '{alias}' via server {server_url}: {e}"
    else:
        return f"Could not resolve genome alias '{alias}' via server {server_url}."
    try:
        listed = client.alias.list(params={"name": alias})
    except Exception as e:  # noqa: BLE001 - the listing only refines the message
        logger.debug(f"Could not list alias '{alias}' on {server_url}: {e}")
        listed = []
    digest = next((a.get("genome_digest") for a in listed if a.get("name") == alias), None)
    if not digest:
        return f"Server {server_url} does not know genome alias '{alias}'."
    return (
        f"Server {server_url} lists genome '{alias}' ({digest}) but did not return "
        f"its sequence collection, so the genome cannot be set up here. This is a "
        f"problem on the server, which is probably not loading the refget store that "
        f"holds this genome. Please report it to the server's maintainers."
    )


_KEPT_EXISTING = "Skipped '{}': existing content was kept. Use --force to replace it."


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


@dataclass(frozen=True)
class PullOptions:
    """
    How a pull treats existing content and large archives.

    Attributes:
        force: Existing content: True replaces it, False keeps it, None asks.
        force_large: Archives over ``size_cutoff``: True pulls, False skips, None asks.
        size_cutoff: The size cutoff in GB; None means ``DEFAULT_PULL_SIZE_CUTOFF_GB``.
        confirm: Answers the questions. Defaults to a refusal unless the CLI
            has enabled interactive prompts; see ``refgenie.utils.prompt``.
    """

    force: bool | None = None
    force_large: bool | None = None
    size_cutoff: int | float | None = None
    confirm: Confirmer | None = None


@dataclass
class BulkPullResult:
    """What a bulk pull did with each asset entry."""

    pulled: list["Asset"] = field(default_factory=list)
    skipped: list[dict] = field(default_factory=list)
    failed: list[dict] = field(default_factory=list)


@dataclass
class PulledAssetMetadata:
    """Server metadata resolved for a single asset pull."""

    asset_digest: str
    declared_modes: list[str]
    asset_parents: list
    asset_name: str
    asset_class_name: str
    asset_names: list
    #: Colocation entries, each with a ``parent_asset_group``.
    colocate: list[dict[str, str]] | None = None
    #: Seek key values from the server, in ``custom_seek_keys`` form.
    seek_keys: dict[str, tuple[str, SeekKeyType]] = field(default_factory=dict)
    #: The server's id of the asset group's class, if it sent one.
    asset_class_id: int | None = None


@dataclass
class PullTransaction:
    """Context manager for atomic pull operations. Rolls back on failure."""

    alias_manager: "AliasManager"
    genome_manager: "GenomeManager"
    _created_genomes: list = field(default_factory=list)
    _created_aliases: list = field(default_factory=list)
    _created_dirs: list = field(default_factory=list)
    #: The genome each created alias names, so ``release_genome`` can find them.
    _alias_genomes: dict = field(default_factory=dict)

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
        """Record a genome/alias creation for rollback. May be called once per alias."""
        if result.created_alias and result.alias_name:
            self._created_aliases.append(result.alias_name)
            self._alias_genomes[result.alias_name] = result.genome_digest
        if result.created_genome and result.genome_digest:
            self._created_genomes.append(result.genome_digest)

    def release_genome(self, genome_digest: str) -> None:
        """Stop rolling back a genome (and its aliases) once a parent asset was pulled into it."""
        self._created_genomes = [g for g in self._created_genomes if g != genome_digest]
        self._created_aliases = [
            a for a in self._created_aliases if self._alias_genomes.get(a) != genome_digest
        ]

    def track_directory(self, path: Path):
        """Track a created directory for potential rollback."""
        self._created_dirs.append(path)


class AssetPuller(ResourceManager):
    """Downloads assets from subscribed servers and records them through ``rgc.asset``."""

    def __init__(
        self,
        database_engine: Engine,
        asset: "AssetManager",
        servers: "ServerManager",
        alias_manager: "AliasManager",
        genome_manager: "GenomeManager",
        events: EventSink | None = None,
        asset_class: "AssetClassManager | None" = None,
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
            asset_class: Where an asset class the server sends is registered
                when the client does not have it yet. Defaults to one over
                ``database_engine``.
        """
        super().__init__(database_engine)
        if asset_class is None:
            from refgenie.managers.asset_class import AssetClassManager

            asset_class = AssetClassManager(database_engine)
        self._asset_class = asset_class
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
        Whether the alias tree for an asset is absent for this alias.

        The tree is derived from the catalog, so a missing tree is safe to re-render.

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
        self, asset_list: list[dict], options: PullOptions = PullOptions()
    ) -> BulkPullResult:
        """
        Pull multiple assets, continuing past failures and logging a summary.

        Args:
            asset_list: Asset dicts from ServerManager.list_assets_for_genome().
                        Each dict should have: genome_digest, asset_group_name, asset_name.
            options: How to treat existing content and large archives.

        Returns:
            BulkPullResult: The pulled assets, and the skipped and failed entries.
                An asset already present counts as skipped.
        """
        result = BulkPullResult()

        for asset_info in asset_list:
            genome_digest = asset_info.get("genome_digest")
            asset_group_name = asset_info.get("asset_group_name")
            asset_name = asset_info.get("asset_name")
            server_url = asset_info.get("server_url")

            if not all([genome_digest, asset_group_name]):
                logger.warning(f"Skipping invalid asset entry: {asset_info}")
                result.failed.append(asset_info)
                continue

            logger.info(f"Pulling {genome_digest}/{asset_group_name}:{asset_name or 'default'}")
            try:
                pulled = self.pull(
                    asset_group_name=asset_group_name,
                    genome=GenomeDigest(genome_digest),
                    asset_name=asset_name,
                    options=options,
                    server_urls=[server_url] if server_url else None,
                )
                if pulled is not None:
                    result.pulled.append(pulled)
                else:
                    result.skipped.append(asset_info)
            except AssetExistsError:
                # Often a parent that an earlier entry in this batch pulled.
                logger.info(
                    f"Already present: {genome_digest}/{asset_group_name}:{asset_name or 'default'}"
                )
                result.skipped.append(asset_info)
            except Exception as e:
                logger.error(f"Failed to pull {genome_digest}/{asset_group_name}:{asset_name}: {e}")
                result.failed.append(asset_info)

        summary = (
            f"Pull summary: {len(result.pulled)} succeeded, {len(result.skipped)} skipped, "
            f"{len(result.failed)} failed"
        )
        if result.failed:
            logger.warning(summary)
        else:
            logger.info(summary)

        return result

    # === Pull Methods ===

    def pull(
        self,
        asset_group_name: str,
        genome: GenomeAlias | GenomeDigest,
        asset_name: str | None = None,
        options: PullOptions = PullOptions(),
        server_urls: list[str] | None = None,
        *,
        required_by: tuple[str, ...] = (),
    ):
        """
        Download and register one asset, firing ``pre_pull`` and ``post_pull``.

        ``required_by`` is internal: digests of the children waiting on this
        pull, used to detect parent cycles. ``post_pull`` fires on every
        outcome. See :meth:`_pull` for the arguments.

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
                options=options,
                server_urls=server_urls,
                required_by=required_by,
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
        options: PullOptions = PullOptions(),
        server_urls: list[str] | None = None,
        required_by: tuple[str, ...] = (),
    ):
        """
        Download and register one asset from the first subscribed server that
        can serve it.

        Args:
            asset_group_name: Name of a group of assets to fetch.
            genome: The genome, as a ``GenomeAlias`` or a ``GenomeDigest``. It
                need not exist locally; see :meth:`_resolve_pull_genome`.
            asset_name: Name of particular asset to fetch.
            options: How to treat existing content and large archives.
            server_urls: Servers to try, in order. Defaults to the subscriptions.

        Returns:
            Asset or None: The added asset, or None if the pull was skipped
            (existing content was kept, or a large archive was declined).

        Raises:
            NoSubscriptionsError: There is no server to try.
            PullFailedError: Every server was tried and none could serve the
                asset. ``reasons`` holds each server's answer. (A
                ``NoArchiveError`` when a server said it has no archive.)
        """
        added_asset = None
        last_error = None
        # Why each server could not serve the asset, reported on failure.
        reasons: dict[str, str] = {}

        server_urls = server_urls or self._servers.subscriptions()
        if not server_urls:
            raise NoSubscriptionsError()

        for server_url in server_urls:
            try:
                added_asset = self._pull_from_server(
                    server_url=server_url,
                    asset_group_name=asset_group_name,
                    genome=genome,
                    asset_name=asset_name,
                    options=options,
                    required_by=required_by,
                )
            except ServerCannotServe as exc:
                logger.warning(f"Server {server_url} cannot serve this asset: {exc}")
                reasons[server_url] = str(exc)
                if exc.final_error is not None:
                    last_error = exc.final_error
                continue
            except PullSkipped as exc:
                logger.warning(str(exc))
                return None

            break

        if added_asset is None:
            if last_error is not None:
                if isinstance(last_error, PullFailedError):
                    last_error.reasons = reasons
                raise last_error
            target = f"{genome}/{asset_group_name}" + (f":{asset_name}" if asset_name else "")
            raise PullFailedError(
                f"Could not pull '{target}': no server could serve it. "
                + " ".join(f"[{url}] {reason}" for url, reason in reasons.items()),
                reasons=reasons,
            )
        return added_asset

    def _pull_from_server(
        self,
        server_url: str,
        asset_group_name: str,
        genome: GenomeAlias | GenomeDigest,
        asset_name: str | None,
        options: PullOptions,
        required_by: tuple[str, ...] = (),
    ) -> "Asset | None":
        """
        One pull attempt against one server, inside its own transaction.

        Raises:
            ServerCannotServe: This server cannot serve the asset; the
                transaction has rolled back and the caller may try the next one.
            PullSkipped: The user or ``force=False`` declined to replace
                existing content, or a large archive was declined.
            MissingAssetClassError: The asset's class is not registered here
                and this server could not provide it.
            ParentPullError: A parent the asset needs is not here and could
                not be pulled; the message names the parent and why.
        """
        with PullTransaction(self._alias_manager, self._genome_manager) as txn:
            client = self._servers.client(server_url)
            alias_name, genome_digest, genome_attrs = self._resolve_pull_genome(
                genome, client, server_url, txn
            )
            try:
                asset_name = str(
                    asset_name
                    or self._asset.group.get_default(asset_group_name, genome_digest=genome_digest)
                )
            except MissingAssetGroupError:
                pass
            else:
                if not options.force and self._asset.exists(
                    genome_digest=genome_digest,
                    asset_group_name=asset_group_name,
                    asset_name=asset_name,
                ):
                    # A pull killed between commit and tree rendering.
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
            # Adopt the server's canonical asset name.
            asset_name = meta.asset_name
            asset_digest = meta.asset_digest

            bundle_name = f"{alias_name}/{asset_group_name}:{asset_name}"

            # The asset class and every parent must be here before downloading.
            meta.asset_class_name = self._ensure_asset_class(client, server_url, meta)
            self._pull_missing_parents(
                parents=meta.asset_parents,
                server_url=server_url,
                child_label=f"{asset_group_name}:{asset_name}",
                child_digest=asset_digest,
                required_by=required_by,
                txn=txn,
                options=options,
            )

            # Content is addressed by digest; the archive's top-level entry is that digest.
            asset_dir = self._asset.data_folder / genome_digest / asset_group_name / asset_digest
            self._clear_existing_content(asset_dir, asset_digest, bundle_name, options)

            # After directory cleanup, so the ORM delete handler cannot remove
            # newly downloaded files.
            if options.force and self._asset.exists(
                genome_digest=genome_digest,
                asset_group_name=asset_group_name,
                asset_name=asset_name,
            ):
                self._asset.remove(
                    genome_digest=genome_digest,
                    asset_group_name=asset_group_name,
                    asset_name=asset_name,
                )

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
                options=options,
            )

            # Inside the transaction so rollback covers it.
            return self._finalize_pulled_asset(
                alias_name=alias_name,
                genome_digest=genome_digest,
                asset_group_name=asset_group_name,
                asset_name=asset_name,
                asset_dir=asset_dir,
                meta=meta,
                genome_attrs=genome_attrs,
            )

    def _ensure_asset_class(
        self, client: "RefgenieserverClient", server_url: str, meta: PulledAssetMetadata
    ) -> str:
        """
        The local name of the asset's class, registering it from the server if needed.

        A fresh install has no asset classes, so the definition is taken from
        the server via the asset group's ``asset_class_id``.

        Returns:
            str: The asset class name to register the pulled asset under.

        Raises:
            MissingAssetClassError: The class is not registered here and the
                server could not provide it.
        """
        name = meta.asset_class_name
        if self._asset_class.exists_any_version(name):
            return name
        if meta.asset_class_id is None:
            raise MissingAssetClassError(name)
        try:
            remote = client.asset_class.get(meta.asset_class_id)
        except Exception as e:
            logger.warning(f"Could not fetch asset class '{name}' from {server_url}: {e}")
            raise MissingAssetClassError(name) from e
        name = remote["name"]
        if self._asset_class.exists_any_version(name):
            return name
        seek_keys = remote.get("seek_keys")
        definition = {
            "name": name,
            "version": remote["version"],
            "description": remote.get("description"),
            "serving_modes": remote.get("serving_modes") or ["file"],
            "seek_keys": {
                sk["name"]: {
                    "value": sk.get("value"),
                    "description": sk.get("description"),
                    "type": sk["type"],
                }
                for sk in seek_keys or []
            },
        }
        self._asset_class.add_from_dict(definition)
        logger.info(f"Registered asset class '{name}' v{remote['version']} from {server_url}")
        if seek_keys is None:
            # Older servers send no class seek keys; only local builds are affected.
            logger.warning(
                f"{server_url} did not send the seek keys of asset class '{name}'. "
                f"Pulled assets are unaffected; to build '{name}' locally, sync the "
                f"class from a data channel with --exists-overwrite."
            )
        return name

    def _pull_missing_parents(
        self,
        parents: list[dict],
        server_url: str,
        child_label: str,
        child_digest: str,
        required_by: tuple[str, ...],
        txn: PullTransaction,
        options: PullOptions,
    ) -> None:
        """
        Pull, from the same server, every parent of an asset that is not here yet.

        Each parent is a full pull with its own transaction, so a pulled parent
        stays even if the child then fails.

        Raises:
            ParentPullError: A parent could not be pulled. The message names it
                (genome, asset group, asset name, digest) and gives the reason.
        """
        chain = (*required_by, child_digest)
        for parent in parents or []:
            digest = parent["digest"]
            if self._asset.exists_digest(digest):
                continue
            genome_digest = parent.get("genome_digest")
            group = parent.get("asset_group_name")
            name = parent.get("name")
            label = f"{genome_digest}/{group}:{name} (digest {digest})"
            if digest in chain:
                raise ParentPullError(
                    f"Cannot pull parent {label} required by {child_label}: the server "
                    f"reports a parent cycle."
                )
            if not genome_digest or not group:
                raise ParentPullError(
                    f"Cannot pull parent {label} required by {child_label}: the server "
                    f"did not say which genome and asset group it belongs to."
                )
            logger.info(
                f"Pulling parent {group}:{name} (digest {digest}) required by {child_label}"
            )
            try:
                pulled = self.pull(
                    asset_group_name=group,
                    genome=GenomeDigest(genome_digest),
                    asset_name=name,
                    # Never replace local content to make room for a parent.
                    options=replace(options, force=False),
                    server_urls=[server_url],
                    required_by=chain,
                )
            except ParentPullError:
                raise
            except Exception as e:
                raise ParentPullError(
                    f"Could not pull parent {label} required by {child_label}: {e}"
                ) from e
            if pulled is None:
                raise ParentPullError(
                    f"Could not pull parent {label} required by {child_label}: the pull "
                    f"was skipped."
                )
            if not self._asset.exists_digest(digest):
                raise ParentPullError(
                    f"Could not pull parent {label} required by {child_label}: the pulled "
                    f"asset has digest {pulled.digest}, not {digest}."
                )
            txn.release_genome(genome_digest)

    def _clear_existing_content(
        self,
        asset_dir: Path,
        asset_digest: str,
        bundle_name: str,
        options: PullOptions,
    ) -> None:
        """
        Make way for a download into ``asset_dir``, or raise ``PullSkipped``.

        Orphaned content (no complete catalog row) is always reclaimed; real
        content is replaced only when ``options.force`` or the confirmer allows it.
        """
        if asset_dir.exists():
            logger.info(f"Asset directory already exists: {asset_dir}")
            try:
                existing_row = self._asset.get_by_digest(asset_digest)
            except MissingAssetError:
                existing_row = None
            if existing_row is None or existing_row.path is None:
                # Debris from a pull killed between extraction and commit.
                logger.warning(
                    f"Reclaiming an orphaned asset directory with no catalog entry: {asset_dir}"
                )
                shutil.rmtree(asset_dir)
            elif options.force is False:
                raise PullSkipped(_KEPT_EXISTING.format(bundle_name))
            elif options.force is None:
                # Never a bare prompt: server and worker pulls have no terminal.
                if not resolve_confirmer(options.confirm)(f"Replace existing {asset_dir}?"):
                    raise PullSkipped(_KEPT_EXISTING.format(bundle_name))
                else:
                    logger.debug(f"Overwriting: {asset_dir}")
                    shutil.rmtree(asset_dir)
            else:
                logger.debug(f"Overwriting: {asset_dir}")
                shutil.rmtree(asset_dir)

    def _download_staged(
        self,
        client: "RefgenieserverClient",
        server_url: str,
        asset_digest: str,
        declared_modes: list[str],
        genome_digest: GenomeDigest,
        asset_group_name: str,
        asset_name: str,
        asset_dir: Path,
        bundle_name: str,
        txn: PullTransaction,
        options: PullOptions,
    ) -> None:
        """
        Download whatever the server has staged for ``asset_digest`` into ``asset_dir``.

        Raises:
            ServerCannotServe: Nothing downloadable is staged on this server.
            PullFailedError: The asset is metadata-only, or staging is ambiguous.
        """
        try:
            staged_items = client.staged_asset.list(params={"asset_digest": asset_digest})
        except Exception as e:
            raise ServerCannotServe(f"Failed to query staged assets from {server_url}: {e}") from e

        archive_records = [rec for rec in staged_items if rec["mode"] == "archive"]
        file_records = [rec for rec in staged_items if rec["mode"] == "file"]

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
                options=options,
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
        client: "RefgenieserverClient",
        server_url: str,
        txn: PullTransaction,
    ) -> tuple[GenomeAlias, GenomeDigest, dict]:
        """
        The alias, the digest and the server's genome record for a pull of ``genome``.

        The alias tree is keyed by alias, so a pull needs both. An unknown
        ``GenomeAlias`` is resolved by the server and registered locally; a
        ``GenomeDigest`` is registered under the server's aliases for it.
        Registrations are tracked on ``txn`` for rollback, and the server's
        genome record fills any unset local genome fields.

        Raises:
            ServerCannotServe: If this server cannot identify the genome.
            PullFailedError: If the server's alias listing fails.
        """
        if isinstance(genome, GenomeAlias):
            alias_name, genome_digest, attrs = self._resolve_pull_alias(
                genome, client, server_url, txn
            )
        else:
            alias_name, genome_digest, attrs = self._resolve_pull_digest(
                genome, client, server_url, txn
            )
        if attrs is None:
            attrs = self._servers.genome_attrs(genome_digest, [server_url])
        self._genome_manager.fill_from_server(genome_digest, attrs)
        return alias_name, genome_digest, attrs

    def _resolve_pull_alias(
        self,
        genome_alias: GenomeAlias,
        client: "RefgenieserverClient",
        server_url: str,
        txn: PullTransaction,
    ) -> tuple[GenomeAlias, GenomeDigest, dict | None]:
        """The pull target for an alias: a local one, else one the server knows.

        The third item is the server's genome record, if it was fetched here.
        """
        try:
            return genome_alias, self._alias_manager.resolve(genome_alias), None
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
            raise ServerCannotServe(_unresolved_alias_reason(genome_alias, client, server_url))
        txn.track_creation(result)
        return genome_alias, self._alias_manager.resolve(genome_alias), result.genome_attrs

    def _resolve_pull_digest(
        self,
        genome_digest: GenomeDigest,
        client: "RefgenieserverClient",
        server_url: str,
        txn: PullTransaction,
    ) -> tuple[GenomeAlias, GenomeDigest, dict]:
        """The pull target for a digest: the server's aliases for it, else a local one.

        The third item is the server's genome record.
        """
        known_digest = self._genome_manager.exists(genome_digest)
        try:
            aliases_metadata = client.alias.list(params={"genome_digest": genome_digest})
        except Exception as e:
            if not known_digest:
                raise ServerCannotServe(
                    f"Could not resolve genome digest '{genome_digest}' via server {server_url}."
                ) from e
            raise PullFailedError(f"Failed to fetch aliases metadata from {server_url}: {e}") from e
        attrs = self._servers.genome_attrs(genome_digest, [server_url])
        for alias_metadata in aliases_metadata:
            result = self._genome_manager.ensure_from_remote(
                alias_name=GenomeAlias(str(alias_metadata["name"])),
                genome_digest=GenomeDigest(alias_metadata["genome_digest"]),
                genome_description=alias_metadata.get("genome_description"),
                genome_attrs=attrs,
            )
            txn.track_creation(result)
        if aliases_metadata:
            return GenomeAlias(str(aliases_metadata[0]["name"])), genome_digest, attrs
        local_aliases = self._alias_manager.get_for_genome(genome_digest) if known_digest else []
        if local_aliases:
            return GenomeAlias(local_aliases[0]), genome_digest, attrs
        raise ServerCannotServe(
            f"Could not resolve genome digest '{genome_digest}' via server {server_url}."
        )

    def _fetch_asset_metadata(
        self,
        client: "RefgenieserverClient",
        server_url: str,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        asset_name: str | None,
        alias_name: GenomeAlias | None,
    ) -> PulledAssetMetadata:
        """
        Query the server for asset group, asset, and relationship metadata,
        including the canonical asset name.
        """
        logger.info(
            f"Querying server {server_url} for {genome_digest}/{asset_group_name}"
            f"{':' + asset_name if asset_name else ''}"
        )

        asset_groups_items = _ask(
            server_url,
            "asset groups",
            client.asset_group.list,
            params={"genome_digest": genome_digest, "asset_group_name": asset_group_name},
        )

        if len(asset_groups_items) != 1:
            raise PullFailedError(
                f"Expected one asset group, got {len(asset_groups_items)} for "
                f"{genome_digest}/{asset_group_name} on server {server_url}"
            )
        asset_group_metadata = asset_groups_items[0]
        asset_group_id = asset_group_metadata["id"]

        assets_items = _ask(
            server_url,
            "assets",
            client.asset.list,
            params={"asset_group_id": asset_group_id, "name": asset_name},
        )

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
        declared_modes = asset_metadata["serving_modes"]

        relationship_metadata = _ask(
            server_url, "asset relationships", client.asset.relationships, asset_digest
        )
        asset_parents = relationship_metadata["parents"]
        asset_name = str(asset_metadata["name"])
        # Older servers omit asset_class_name; fall back to the group name.
        asset_class_name = asset_metadata.get("asset_class_name") or asset_group_metadata["name"]
        asset_names = asset_metadata["names"]

        return PulledAssetMetadata(
            asset_digest=asset_digest,
            declared_modes=declared_modes,
            asset_parents=asset_parents,
            asset_name=asset_name,
            asset_class_name=asset_class_name,
            asset_names=asset_names,
            colocate=asset_metadata.get("colocate"),
            seek_keys=_seek_keys_from_server(asset_metadata.get("seek_keys")),
            asset_class_id=asset_group_metadata.get("asset_class_id"),
        )

    def _finalize_pulled_asset(
        self,
        alias_name: GenomeAlias,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        asset_name: str,
        asset_dir: Path,
        meta: PulledAssetMetadata,
        genome_attrs: dict,
    ):
        """
        Register a downloaded asset: colocation symlinks, fasta genome init,
        catalog write, extra names, default asset, and parent links.

        Must be called inside the PullTransaction.
        """
        progress.emit("stage", message="Registering asset", phase="register")
        asset_digest = meta.asset_digest
        asset_class_name = meta.asset_class_name
        asset_names = meta.asset_names
        asset_parents = meta.asset_parents

        colocate_metadata = meta.colocate
        if colocate_metadata:
            recreate_colocation_symlinks(
                output_folder=asset_dir,
                genome_folder=self._asset.genome_folder,
                colocate_metadata=colocate_metadata,
                parent_assets=self._resolve_colocation_parents(colocate_metadata, genome_digest),
            )

        if asset_class_name == "fasta":
            # expected_digest rejects a bad fasta before any genome row is written.
            try:
                self._genome_manager.initialize_genome(
                    fasta_file_path=asset_dir / f"{genome_digest}.fa",
                    description=genome_attrs.get("description") or "",
                    alias_names=[alias_name],
                    use_existing=True,
                    expected_digest=genome_digest,
                )
            except FastaDigestMismatchError as e:
                raise PullFailedError(
                    f"The fasta the server sent for '{alias_name}' does not match "
                    f"its genome ({e}). The server's fasta asset is probably "
                    f"corrupt; please report it. If your refgenie is old, "
                    f"`pip install -U refgenie` may also fix it."
                ) from e
            # The genome may be new, and the fasta brings what the FHR sidecar needs.
            self._genome_manager.fill_from_server(genome_digest, genome_attrs)

        # Use the server's digest: the recipe's `inherent` set, which decides
        # what the digest covers, does not travel with the archive.
        progress.emit("stage", message="Registering and linking", phase="symlink")
        added_asset = self._asset.content.add(
            asset_class_name=asset_class_name,
            path=asset_dir.relative_to(self._asset.genome_folder),
            genome_digest=genome_digest,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
            colocate=colocate_metadata,
            digest=asset_digest,
            custom_seek_keys=meta.seek_keys,
            build_provenance=_provenance_from_server(asset_names, asset_name),
        )

        # Adopt the server's other names; a conflicting one is skipped with a warning.
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

        try:
            default_asset = self._asset.group.get_default(
                asset_group_name, genome_digest=genome_digest
            )
        except MissingAssetGroupError:
            default_asset = None
        if default_asset is None:
            self._asset.group.set_default(asset_group_name, asset_name, genome_digest=genome_digest)

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
        client: "RefgenieserverClient",
        staged_asset_metadata: dict,
        asset_digest: str,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        asset_name: str,
        asset_dir: Path,
        bundle_name: str,
        txn: PullTransaction,
        options: PullOptions,
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
            options: How to treat large archives.
        """
        if not should_pull_large_archive(
            archive_size=staged_asset_metadata.get("tarball_size"),
            asset_registry_path=bundle_name,
            force=options.force_large,
            size_cutoff=options.size_cutoff,
            confirm=options.confirm,
        ):
            raise PullSkipped(
                f"Skipped '{bundle_name}': its archive is over the size cutoff. "
                f"Use --pull-large or a higher --size-cutoff to pull it."
            )

        # The archive's top-level entry is asset_dir's basename, so download
        # and untar into its parent.
        tardir = asset_dir.parent
        tarpath = tardir / f"{asset_group_name}__{asset_name}.tgz"

        tardir_existed = tardir.exists()
        if not tardir_existed:
            logger.debug(f"Creating directory: {tardir}")
            tardir.mkdir(parents=True, exist_ok=True)
            txn.track_directory(tardir)

        # The download removes its own ``.part`` file on the guard's SystemExit.
        with sigint_guard(f"The download was interrupted: {tarpath}"):
            client.download_with_progress(
                operation_id=API_ID_ARCHIVE,
                output_path=tarpath,
                url_format_params={"asset_digest": asset_digest},
                name=bundle_name,
            )
        # `tarball_digest` is the digest of the archive bytes, not the asset digest.
        # Hashing a large archive is slow, so announce it.
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
        progress.emit("stage", message=f"Unpacking {bundle_name}", phase="unpack")
        logger.info(f"Extracting asset tarball: {tarpath}")
        untar(tarpath.as_posix(), tardir.as_posix(), filter="fully_trusted")
        tarpath.unlink()

    def _pull_file_mode(
        self,
        client: "RefgenieserverClient",
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
        file_list = client.asset.files(asset_digest)
        if not file_list:
            msg = f"No files found for asset {asset_digest}"
            logger.error(msg)
            raise RuntimeError(msg)

        logger.info(f"Downloading {len(file_list)} file(s) for '{bundle_name}'")

        if not asset_dir.exists():
            asset_dir.mkdir(parents=True, exist_ok=True)
            txn.track_directory(asset_dir)

        for file_path in file_list:
            output_path = asset_dir / file_path
            output_path.parent.mkdir(parents=True, exist_ok=True)
            client.asset.download_file(asset_digest, file_path, output_path)

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
