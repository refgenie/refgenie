"""Pulling from subscribed servers: one asset, many genomes, or a full mirror."""

from typing import TYPE_CHECKING

from refgenie.const import DEFAULT_SERVER_URL
from refgenie.exceptions import MissingAliasError, NoSubscriptionsError, PullFailedError
from refgenie.logger import logger
from refgenie.managers.sources.servers import estimate_pull_size
from refgenie.managers.transfer.puller import BulkPullResult, PullOptions
from refgenie.models import GenomeAlias, GenomeDigest
from refgenie.plugins.events import EventSink, update_scope
from refgenie.utils.prompt import Confirmer, resolve_confirmer

if TYPE_CHECKING:
    from refgenie.db.tables import Asset
    from refgenie.managers.alias import AliasBackend
    from refgenie.managers.genome import GenomeManager
    from refgenie.managers.sources.servers import ServerManager
    from refgenie.managers.transfer.puller import AssetPuller


def confirm_bulk_pull(
    asset_count: int,
    genome_count: int,
    total_bytes: int,
    force: bool = False,
    confirm: Confirmer | None = None,
) -> bool:
    """
    Display size estimate and ask for confirmation.

    Args:
        asset_count: Number of assets to pull.
        genome_count: Number of genomes involved.
        total_bytes: Total estimated download size in bytes.
        force: If True, return True without prompting.
        confirm: Confirmation callback. Defaults to a refusal unless the CLI
            has enabled interactive prompts; see `refgenie.utils.prompt`.

    Returns:
        True if user confirms or force=True, False otherwise.
    """
    if force:
        return True

    if total_bytes > 0:
        size_gb = total_bytes / (1024**3)
        size_str = f"{size_gb:.1f} GB"
    else:
        size_str = "unknown size"

    msg = (
        f"About to pull {asset_count} asset(s) for {genome_count} genome(s) "
        f"(estimated {size_str}). Continue?"
    )
    return resolve_confirmer(confirm)(msg)


class TransferManager:
    """
    Pulling assets from subscribed servers: one asset, many genomes, or a full mirror.

    Each public method is one ``update_scope``, so a bulk pull of N assets fires
    N pre/post pull pairs and one ``post_update``.
    """

    def __init__(
        self,
        servers: "ServerManager",
        alias_manager: "AliasBackend",
        genome_manager: "GenomeManager",
        puller: "AssetPuller",
        events: EventSink,
    ):
        self._servers = servers
        self._alias = alias_manager
        self._genome = genome_manager
        self._puller = puller
        self._events = events  # required by @update_scope

    @update_scope
    def pull(
        self,
        asset_group_name: str,
        genome: GenomeAlias | GenomeDigest,
        asset_name: str | None = None,
        force: bool | None = None,
        force_large: bool | None = None,
        force_server_urls: list[str] | None = None,
        size_cutoff: int | float | None = None,
        confirm: Confirmer | None = None,
    ) -> "Asset | None":
        """
        Download and unpack an asset for a given reference genome.

        Args:
            asset_group_name: Name of a group of assets to fetch.
            genome: The genome, as a ``GenomeAlias`` or a ``GenomeDigest``. It
                need not exist locally; the servers are asked about it. A
                plain ``str`` raises TypeError: say which one it is.
            asset_name: Name of particular asset to fetch.
            force: How to handle case in which asset path already exists.
            force_large: How to handle archives larger than size_cutoff (default 10GB).
            force_server_urls: Force specific server URLs to use.
            size_cutoff: Maximum archive file size to download without prompt.
            confirm: Confirmation callback. Defaults to a refusal unless the CLI
                has enabled interactive prompts; see `refgenie.utils.prompt`.

        Returns:
            Asset or None: The added asset, or None if the pull was skipped
            (existing content was kept, or a large archive was declined).

        Raises:
            NoSubscriptionsError: Nothing is subscribed and the offer to
                subscribe to the default server was declined.
            PullFailedError: No server could serve the asset; ``reasons``
                says why, per server.
        """
        if not force_server_urls and not self._servers.subscriptions():
            logger.error("No server subscriptions found")
            if not resolve_confirmer(confirm)("Would you like to subscribe to the default server?"):
                raise NoSubscriptionsError()
            self._servers.subscribe([DEFAULT_SERVER_URL])

        return self._puller.pull(
            asset_group_name=asset_group_name,
            genome=genome,
            asset_name=asset_name,
            options=PullOptions(force, force_large, size_cutoff, confirm),
            server_urls=force_server_urls,
        )

    @update_scope
    def pull_genomes(
        self,
        genomes: list[GenomeAlias | GenomeDigest] | None = None,
        asset_group_name: str | None = None,
        all_genomes: bool = False,
        force: bool | None = None,
        force_large: bool | None = None,
        size_cutoff: int | float | None = None,
        confirm: Confirmer | None = None,
    ) -> BulkPullResult:
        """
        Pull every asset, or every asset in one asset group, for many genomes.

        Shows a size estimate and confirmation prompt (unless force=True).

        Args:
            genomes: Genomes, each a ``GenomeAlias`` or a ``GenomeDigest``.
                Ignored if all_genomes=True.
            asset_group_name: Only pull assets in this asset group (e.g. "fasta").
                None pulls every asset.
            all_genomes: If True, queries servers for all available genomes.
            force: Skip confirmation prompts if True.
            force_large: How to handle large archives.
            size_cutoff: Maximum archive file size to download without prompt.
            confirm: Confirmation callback.

        Returns:
            BulkPullResult: What was pulled, skipped and failed.
        """
        if all_genomes:
            genomes = self._remote_genome_refs()

        if not genomes:
            logger.error("No genomes specified")
            return BulkPullResult()

        digests = self._resolve_or_init(genomes)
        if not digests:
            logger.error("No valid genomes to pull")
            return BulkPullResult()

        empty_message = (
            "No remote assets found for specified genomes"
            if asset_group_name is None
            else f"No remote assets named '{asset_group_name}' found for specified genomes"
        )
        return self._confirm_and_pull(
            self._remote_assets(digests, asset_group_name),
            len(digests),
            PullOptions(force, force_large, size_cutoff, confirm),
            empty_message=empty_message,
            cancelled_message="Bulk pull cancelled",
        )

    @update_scope
    def init_genomes(
        self,
        aliases: list[GenomeAlias] | None = None,
        all_genomes: bool = False,
    ) -> list[GenomeDigest]:
        """
        Register genome(s) in the local store from remote metadata, without downloading files.

        Args:
            aliases: Genome aliases to register.
            all_genomes: If True, register all genomes from remote servers
                (genomes with no alias on the server are skipped).

        Returns:
            List of registered genome digests.

        Raises:
            PullFailedError: Some genomes could not be registered (the rest
                still are); the message names them.
        """
        if all_genomes:
            aliases = [
                GenomeAlias(g["aliases"][0])
                for g in self._servers.list_genomes()
                if g.get("aliases")
            ]
        registered, failed = [], []
        for alias in aliases or []:
            if not self._genome.init_from_remote(alias, require_collection=True):
                failed.append(alias)
                continue
            try:
                registered.append(self._alias.resolve(alias))
            except MissingAliasError:
                pass

        logger.info(f"Registered {len(registered)} genome(s)")
        if failed:
            raise PullFailedError(f"could not register {', '.join(failed)}")
        return registered

    @update_scope
    def mirror(
        self,
        force: bool | None = None,
        force_large: bool | None = None,
        size_cutoff: int | float | None = None,
        confirm: Confirmer | None = None,
    ) -> BulkPullResult:
        """
        Mirror all assets for all genomes from subscribed servers.

        Always shows confirmation with total size estimate.

        Args:
            force: Skip confirmation prompts if True.
            force_large: How to handle large archives.
            size_cutoff: Maximum archive file size to download without prompt.
            confirm: Confirmation callback.

        Returns:
            BulkPullResult: What was pulled, skipped and failed.
        """
        remote_genomes = self._servers.list_genomes()
        if not remote_genomes:
            logger.warning("No remote genomes found")
            return BulkPullResult()

        # Register genomes first. One with no server alias uses its digest as
        # its alias, so the pull has an alias tree to render.
        for g in remote_genomes:
            alias_name = g["aliases"][0] if g.get("aliases") else None
            self._genome.init_from_remote(
                alias_name=GenomeAlias(alias_name or g["genome_digest"]),
                genome_digest=GenomeDigest(g["genome_digest"]),
                genome_description=g.get("description"),
            )

        return self._confirm_and_pull(
            self._remote_assets([GenomeDigest(g["genome_digest"]) for g in remote_genomes]),
            len(remote_genomes),
            PullOptions(force, force_large, size_cutoff, confirm),
            empty_message="No remote assets found to mirror",
            cancelled_message="Mirror cancelled",
        )

    def _resolve_or_init(self, genomes: list[GenomeAlias | GenomeDigest]) -> list[GenomeDigest]:
        """
        Digests of ``genomes``, registering from the servers any alias not known here.

        An alias no server knows is logged and skipped.
        """
        digests = []
        for genome in genomes:
            if isinstance(genome, GenomeDigest):
                digests.append(genome)
                continue
            if not isinstance(genome, GenomeAlias):
                raise TypeError(
                    f"Expected a GenomeAlias or a GenomeDigest, not {type(genome).__name__} "
                    f"'{genome}'."
                )
            try:
                digests.append(self._alias.resolve(genome))
                continue
            except MissingAliasError:
                pass
            if self._genome.init_from_remote(genome):
                digests.append(self._alias.resolve(genome))
            else:
                logger.error(f"Could not resolve genome: {genome}")
        return digests

    def _remote_genome_refs(self) -> list[GenomeAlias | GenomeDigest]:
        """Every genome the servers list: its first alias, or its digest if it has none."""
        return [
            GenomeAlias(g["aliases"][0]) if g["aliases"] else GenomeDigest(g["genome_digest"])
            for g in self._servers.list_genomes()
        ]

    def _remote_assets(
        self, digests: list[GenomeDigest], asset_group_name: str | None = None
    ) -> list[dict]:
        """Remote assets for ``digests``, only those in ``asset_group_name`` if given."""
        assets = []
        for digest in digests:
            for asset in self._servers.list_assets_for_genome(digest):
                if asset_group_name is None or asset.get("asset_group_name") == asset_group_name:
                    assets.append(asset)
        return assets

    def _confirm_and_pull(
        self,
        assets: list[dict],
        genome_count: int,
        options: PullOptions,
        *,
        empty_message: str,
        cancelled_message: str,
    ) -> BulkPullResult:
        """Show the size estimate, ask once, then pull ``assets``."""
        if not assets:
            logger.warning(empty_message)
            return BulkPullResult()

        total_bytes, asset_count = estimate_pull_size(assets)
        if not confirm_bulk_pull(
            asset_count=asset_count,
            genome_count=genome_count,
            total_bytes=total_bytes,
            force=options.force or False,
            confirm=options.confirm,
        ):
            logger.info(cancelled_message)
            return BulkPullResult()

        return self._puller.pull_multiple(asset_list=assets, options=options)
