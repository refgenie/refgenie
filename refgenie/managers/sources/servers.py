"""
ServerManager - the refgenieservers this node pulls from.

Reached as ``rgc.servers``. It owns the subscription list (stored on the
``Configuration`` row), one client per server URL, read-only queries of a
server's catalog, the server-side lookups of an alias or a collection, and
remote seek. Pulling itself is ``AssetPuller``'s job; it
asks this manager for clients and subscriptions.

"Remote" in refgenie means a push target (``refgenie remote``, the ``Remote``
table). Those are not handled here.
"""

from collections import defaultdict
from collections.abc import Callable
from typing import TYPE_CHECKING

from rich.table import Table
from sqlalchemy.engine import Engine

from refgenie.db.tables import is_path_type
from refgenie.exceptions import (
    MissingAliasError,
    MissingAssetError,
    MissingAssetGroupError,
    MissingSeekKeyError,
)
from refgenie.logger import logger
from refgenie.managers.asset.seek_key import default_seek_key_from_payload, payload_is_path_type
from refgenie.managers.asset.tables import asset_table
from refgenie.managers.base import ResourceManager
from refgenie.managers.configuration import latest_configuration
from refgenie.managers.sources.api_ids import API_ID_ALIAS_DIGEST, API_ID_GENOME_ATTRS
from refgenie.managers.sources.client import RefgenieserverClient, ServerClient
from refgenie.managers.sources.genomes import (
    RemoteGenomeSource,
    make_source,
    normalize_server_url,
)
from refgenie.models import AssetRegistryPathComponents, GenomeAlias, GenomeDigest

if TYPE_CHECKING:
    from refgenie.managers.asset.seek_key import SeekKeyManager
    from refgenie.managers.alias import AliasBackend


def estimate_pull_size(asset_list: list[dict]) -> tuple[int, int]:
    """
    Total size and count for a list of remote assets.

    Args:
        asset_list: Asset dicts from ``ServerManager.list_assets_for_genome``,
            each with an ``archive_size`` key.

    Returns:
        (total_bytes, asset_count). Assets of unknown size are counted but add
        nothing to the total.
    """
    total_bytes = sum(a["archive_size"] for a in asset_list if a.get("archive_size") is not None)
    return total_bytes, len(asset_list)


class ServerManager(ResourceManager):
    """
    Refgenieservers this node pulls from: subscriptions, clients, catalog, remote seek.

    It is also the one place that asks servers what an alias or a collection
    is (``resolve_alias``, ``find_collection``, ``genome_source``); genome
    creation, ``set_genome_alias``, remote seek and ``refgenie id --remote``
    all go through it.
    """

    def __init__(
        self,
        database_engine: Engine,
        alias_manager: "AliasBackend",
        seek_keys: "SeekKeyManager",
        clients: dict[str, ServerClient | None] | None = None,
    ):
        """
        Initialize the ServerManager.

        Args:
            database_engine: The database engine.
            alias_manager: The alias manager, for resolving registry-path aliases.
            seek_keys: The SeekKeyManager, for local non-path seek keys in remote seek.
            clients: Pre-built clients by server URL (tests inject these).
                Each must satisfy the ``ServerClient`` protocol.

        Raises:
            ValueError: If a pre-built client does not satisfy ``ServerClient``.
        """
        super().__init__(database_engine)
        for url, client in (clients or {}).items():
            if not isinstance(client, ServerClient):
                raise ValueError(
                    f"Invalid server client for {url}. Does not match ServerClient Protocol"
                )
        self._alias_manager = alias_manager
        self._seek_keys = seek_keys
        self._clients: dict[str, ServerClient] = dict(clients or {})

    # =========================================================================
    # Subscriptions
    # =========================================================================

    def subscribe(self, server_urls: str | list[str], reset: bool = False) -> None:
        """
        Subscribe to a list of servers.

        URLs are stored in normalized form (see ``normalize_server_url``).

        Args:
            server_urls: Server URL or list of server URLs to subscribe to.
            reset: If True, overwrite the current list instead of adding to it.
        """
        if isinstance(server_urls, str):
            server_urls = [server_urls]
        server_urls = [normalize_server_url(url) for url in server_urls]
        with self._database_session as session:
            configuration = latest_configuration(session)
            configuration.servers = (
                list(server_urls) if reset else list(set(configuration.servers).union(server_urls))
            )
            session.add(configuration)
            session.commit()
        logger.info(f"Subscribed to servers: {' '.join(server_urls)}")

    def unsubscribe(self, server_urls: list[str]) -> None:
        """
        Unsubscribe from a list of servers.

        The given URLs are normalized before comparison, matching how
        ``subscribe`` stores them.

        Args:
            server_urls: The list of server URLs to unsubscribe from.
        """
        if not isinstance(server_urls, list):
            raise TypeError(f"servers must be a list of strings, not {type(server_urls)}")
        removing = {normalize_server_url(url) for url in server_urls}
        with self._database_session as session:
            configuration = latest_configuration(session)
            configuration.servers = [
                stored for stored in configuration.servers if stored not in removing
            ]
            session.add(configuration)
            session.commit()
        logger.info(f"Unsubscribed from servers: {server_urls}")

    def subscriptions(self) -> list[str]:
        """The subscribed server URLs."""
        with self._database_session as session:
            return list(latest_configuration(session).servers)

    def find_subscription(self, url: str) -> str | None:
        """
        The stored subscription that ``url`` names (after normalization), or None.

        Builds no client and does no network I/O: the web layer uses this to
        refuse an unsubscribed URL before anything contacts it.

        Args:
            url: The server URL to look up.

        Returns:
            The stored subscription string, or None if ``url`` is not subscribed.
        """
        wanted = normalize_server_url(url)
        return wanted if wanted in self.subscriptions() else None

    # =========================================================================
    # Clients
    # =========================================================================

    @property
    def clients(self) -> dict[str, ServerClient]:
        """The clients built so far, by server URL."""
        return self._clients

    def client(self, url: str) -> ServerClient:
        """
        The client for a server URL, created on first use and cached.

        Args:
            url: The server URL.

        Returns:
            ServerClient: The server client.
        """
        if url not in self._clients:
            self._clients[url] = RefgenieserverClient(url)
        return self._clients[url]

    # =========================================================================
    # Server-side lookups: what an alias or a collection is
    # =========================================================================

    def resolve_alias(
        self, alias: GenomeAlias, server_urls: list[str] | None = None
    ) -> GenomeDigest | None:
        """
        Ask the servers which genome an alias names. Nothing is written locally.

        Each server is asked through its genome source first
        (``RefgenieServerSource.resolve_alias``), then through the cached v4
        client, for servers that are not usable as a genome source.

        Args:
            alias: The alias to resolve.
            server_urls: Servers to ask, in order. Defaults to the subscriptions.

        Returns:
            The digest from the first server that knows the alias, or None.
        """
        server_urls = server_urls or self.subscriptions()
        for server_url in server_urls:
            try:
                digest = make_source(server_url).resolve_alias(alias)
            except Exception as exc:
                logger.debug(f"Failed to create source for {server_url}: {exc}")
                continue
            if digest is not None:
                logger.info(f"Resolved alias '{alias}' to digest: {digest}")
                return GenomeDigest(digest)
            logger.debug(f"Alias '{alias}' not found on {server_url}")
        for server_url in server_urls:
            try:
                data = self.client(server_url).get(
                    operation_id=API_ID_ALIAS_DIGEST, url_format_params={"name": alias}
                )
            except Exception as exc:
                logger.debug(f"Alias resolution via server client failed for {server_url}: {exc}")
                continue
            digest = data.get("digest") if isinstance(data, dict) else None
            if digest:
                logger.info(f"Resolved alias '{alias}' via server client: {digest}")
                return GenomeDigest(digest)
        return None

    def find_collection(
        self, digest: GenomeDigest, server_urls: list[str] | None = None
    ) -> dict | None:
        """
        Ask the servers for a sequence collection by its digest.

        Args:
            digest: The collection (genome) digest.
            server_urls: Servers to ask, in order. Defaults to the subscriptions.

        Returns:
            The collection's level 2 data from the first server that has it, or None.
        """
        for server_url in server_urls or self.subscriptions():
            try:
                collection = make_source(server_url).verify_collection(digest)
            except Exception as exc:
                logger.debug(f"Collection lookup failed on {server_url}: {exc}")
                continue
            if collection is not None:
                return collection
        return None

    def genome_source(self, server_urls: list[str] | None = None) -> RemoteGenomeSource | None:
        """
        The first server usable as a genome source, to create genomes from.

        Args:
            server_urls: Servers to try, in order. Defaults to the subscriptions.

        Returns:
            The source, or None if no server is usable as one.
        """
        for server_url in server_urls or self.subscriptions():
            try:
                return make_source(server_url)
            except Exception as exc:
                logger.debug(f"Failed to create source for {server_url}: {exc}")
        return None

    def genome_description(
        self, genome_digest: GenomeDigest, server_urls: list[str] | None = None
    ) -> str:
        """
        The first description a server gives for a genome.

        Args:
            genome_digest: The genome to describe.
            server_urls: Servers to ask, in order. Defaults to the subscriptions.

        Returns:
            The description, or an empty string if no server has one.
        """
        for server_url in server_urls or self.subscriptions():
            try:
                attrs = self.client(server_url).get(
                    operation_id=API_ID_GENOME_ATTRS,
                    url_format_params={"digest": genome_digest},
                )
            except Exception as exc:
                logger.debug(f"No genome description from {server_url}: {exc}")
                continue
            if isinstance(attrs, dict) and attrs.get("description"):
                return attrs["description"]
        return ""

    # =========================================================================
    # Catalog: what a server has
    # =========================================================================

    def list_assets(
        self,
        genome_digests: list[GenomeDigest] | None = None,
        include_seek_keys: bool = False,
        server_urls: list[str] | None = None,
    ) -> tuple[dict[str, dict[str, list[str]]], dict[str, dict[str, str]]]:
        """
        List all assets available on the subscribed servers.

        Args:
            genome_digests: Optional list of genome digests to filter by.
            include_seek_keys: Accepted for symmetry with the local listing;
                servers do not report seek keys here, so it changes nothing.
            server_urls: Servers to query. Defaults to the subscriptions.

        Returns:
            (assets by server URL, aliases by server URL). Asset data has the
            shape of ``AssetManager.list_all``: genome digest -> ``group:asset``
            strings. Alias data is genome digest -> comma-separated aliases.
        """
        server_urls = server_urls or self.subscriptions()
        if not server_urls:
            logger.warning("No server subscriptions found")
            return {}, {}

        result = {}
        aliases_result = {}

        for server_url in server_urls:
            logger.debug(f"Querying remote assets from: {server_url}")
            client = self.client(server_url)

            assets_response = client.get_all_assets()
            if not assets_response:
                logger.warning(f"No assets found on server: {server_url}")
                result[server_url] = {}
                aliases_result[server_url] = {}
                continue

            asset_groups_response = client.get_all_asset_groups()
            try:
                aliases_response = client.get_all_aliases()
            except Exception as e:
                logger.warning(f"Failed to fetch aliases from {server_url}: {e}")
                aliases_response = []

            asset_groups = {g["id"]: g for g in asset_groups_response or []}

            server_assets: dict[str, list[str]] = {}
            for asset in assets_response:
                asset_group = asset_groups.get(asset.get("asset_group_id"))
                if not asset_group:
                    continue
                genome_digest = asset_group.get("genome_digest")
                if not genome_digest:
                    continue
                if genome_digests and genome_digest not in genome_digests:
                    continue
                server_assets.setdefault(genome_digest, []).append(
                    f"{asset_group.get('name', 'unknown')}:{asset.get('name', 'unknown')}"
                )

            aliases_by_genome = defaultdict(list)
            for alias_item in aliases_response:
                genome_digest = alias_item.get("genome_digest")
                alias_name = alias_item.get("name")
                if genome_digest in server_assets and alias_name:
                    aliases_by_genome[genome_digest].append(alias_name)

            result[server_url] = server_assets
            aliases_result[server_url] = {
                digest: ", ".join(names) for digest, names in aliases_by_genome.items()
            }
            logger.debug(f"Found {len(server_assets)} genomes with assets on {server_url}")

        return result, aliases_result

    def list_assets_for_genome(
        self,
        genome_digest: GenomeDigest,
        server_urls: list[str] | None = None,
    ) -> list[dict]:
        """
        List all assets available for a genome on the servers.

        Args:
            genome_digest: The digest of the genome to query.
            server_urls: Servers to query. Defaults to the subscriptions.

        Returns:
            Dicts with keys server_url, genome_digest, asset_group_name,
            asset_name, asset_digest, archive_digest, archive_size.
        """
        server_urls = server_urls or self.subscriptions()
        if not server_urls:
            logger.warning("No server subscriptions found")
            return []

        results = []
        for server_url in server_urls:
            logger.debug(f"Querying remote assets for genome {genome_digest} from: {server_url}")
            try:
                client = self.client(server_url)
                asset_groups = client.get_all_asset_groups(params={"genome_digest": genome_digest})
            except Exception as e:
                logger.warning(f"Failed to fetch asset groups from {server_url}: {e}")
                continue

            for asset_group in asset_groups or []:
                asset_group_id = asset_group.get("id")
                asset_group_name = asset_group.get("name", "unknown")

                try:
                    assets = client.get_all_assets(params={"asset_group_id": asset_group_id})
                except Exception as e:
                    logger.warning(f"Failed to fetch assets from {server_url}: {e}")
                    continue

                for asset in assets or []:
                    asset_digest = asset.get("digest")
                    archive_size = None
                    archive_digest = None
                    if asset_digest:
                        try:
                            staged = client.get_staged_assets(params={"asset_digest": asset_digest})
                            if staged:
                                # tarball_* only: a record's `digest` is the asset
                                # identity digest, not the archive's byte digest.
                                archive_size = staged[0].get("tarball_size")
                                archive_digest = staged[0].get("tarball_digest")
                        except Exception as e:
                            logger.debug(f"Failed to fetch archive size for {asset_digest}: {e}")

                    results.append(
                        {
                            "server_url": server_url,
                            "genome_digest": genome_digest,
                            "asset_group_name": asset_group_name,
                            "asset_name": asset.get("name", "unknown"),
                            "asset_digest": asset_digest,
                            "archive_digest": archive_digest,
                            "archive_size": archive_size,
                        }
                    )

        return results

    def list_genomes(self, server_urls: list[str] | None = None) -> list[dict]:
        """
        List all genomes available on the servers, first server wins per digest.

        Args:
            server_urls: Servers to query. Defaults to the subscriptions.

        Returns:
            Dicts with keys server_url, genome_digest, aliases, description.
        """
        server_urls = server_urls or self.subscriptions()
        if not server_urls:
            logger.warning("No server subscriptions found")
            return []

        results = []
        seen_digests = set()

        for server_url in server_urls:
            logger.debug(f"Querying remote genomes from: {server_url}")
            try:
                client = self.client(server_url)
                genomes = client.get_all_genomes()
            except Exception as e:
                logger.warning(f"Failed to fetch genomes from {server_url}: {e}")
                continue

            try:
                aliases = client.get_all_aliases()
            except Exception as e:
                logger.warning(f"Failed to fetch aliases from {server_url}: {e}")
                aliases = []

            aliases_by_genome = defaultdict(list)
            for alias in aliases or []:
                genome_digest = alias.get("genome_digest")
                alias_name = alias.get("name")
                if genome_digest and alias_name:
                    aliases_by_genome[genome_digest].append(alias_name)

            for genome in genomes or []:
                genome_digest = genome.get("digest")
                if genome_digest and genome_digest not in seen_digests:
                    seen_digests.add(genome_digest)
                    results.append(
                        {
                            "server_url": server_url,
                            "genome_digest": genome_digest,
                            "aliases": aliases_by_genome.get(genome_digest, []),
                            "description": genome.get("description", ""),
                        }
                    )

        return results

    def assets_table(
        self,
        genome_digests: list[GenomeDigest] | None = None,
        include_seek_keys: bool = False,
        server_urls: list[str] | None = None,
    ) -> list[Table]:
        """
        Tables of the assets on the servers, one per server.

        Args:
            genome_digests: The digests of genomes to filter by.
            include_seek_keys: Whether to include seek keys.
            server_urls: Servers to query. Defaults to the subscriptions.

        Returns:
            list[Table]: One table per server.
        """
        asset_data_by_source, aliases_data_by_source = self.list_assets(
            genome_digests=genome_digests,
            include_seek_keys=include_seek_keys,
            server_urls=server_urls,
        )
        return [
            asset_table(
                asset_data=asset_data,
                aliases_data=aliases_data_by_source.get(source, {}),
                include_seek_keys=include_seek_keys,
                source=source,
            )
            for source, asset_data in asset_data_by_source.items()
        ]

    # =========================================================================
    # Remote seek: registry path -> file URL on a server
    # =========================================================================

    def seek(
        self,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        asset_name: str | None = None,
        seek_key: str | None = None,
        server_urls: list[str] | None = None,
    ) -> str:
        """
        Seek a remote path to an asset via the file-level endpoint.

        For assets with file serving mode, returns a direct URL to the file.
        For non-path seek keys, returns the value directly.

        Args:
            genome_digest: The genome digest. It need not exist locally.
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            seek_key: The seek key.
            server_urls: Servers to query. Defaults to the subscriptions.

        Returns:
            str: The remote URL to the asset file, or the value of a non-path seek key.
        """
        components = AssetRegistryPathComponents(
            asset_group=asset_group_name, asset=asset_name, seek_key=seek_key
        )
        return self._seek(components, lambda server_url: genome_digest, genome_digest, server_urls)

    def _local_non_path_value(
        self, genome_digest: GenomeDigest, components: AssetRegistryPathComponents
    ) -> str | None:
        """
        The value of a non-path seek key, if this node has the asset.

        A non-path value (a version string, a JSON blob) is the same wherever it
        is read, so a local copy answers without asking a server. Returns None
        when the asset is not here or the seek key is a path.
        """
        try:
            _asset_name, seek_key = self._seek_keys.resolve(
                genome_digest, components.asset_group, components.asset, components.seek_key
            )
        except (MissingAssetGroupError, MissingAssetError, MissingSeekKeyError):
            return None
        return None if is_path_type(seek_key.type) else seek_key.value

    def seek_components(
        self,
        components: AssetRegistryPathComponents,
        server_urls: list[str] | None = None,
    ) -> str:
        """
        Get the remote path to an asset from registry path components.

        The genome in a registry path is an alias. A local alias is used when
        there is one; otherwise each server is asked to resolve it (read-only,
        nothing is written locally). For non-path seek keys, the value is
        returned directly (no remote URL needed). For path-based seek keys,
        constructs a URL to the file-level endpoint.

        Args:
            components: The asset registry path components.
            server_urls: Servers to query. Defaults to the subscriptions.

        Returns:
            str: The remote URL to the asset file, or the value of a non-path seek key.

        Raises:
            ValueError: If the asset does not support file-level access, or no
                server could resolve it.
        """
        if components.genome is None:
            raise ValueError(f"Genome name is required: {components=}")
        genome_alias = GenomeAlias(components.genome)
        try:
            local_digest = self._alias_manager.resolve(genome_alias)
        except MissingAliasError:
            local_digest = None

        def digest_on(server_url: str) -> GenomeDigest | None:
            if local_digest is not None:
                return local_digest
            return self.resolve_alias(genome_alias, [server_url])

        return self._seek(components, digest_on, local_digest, server_urls)

    def _seek(
        self,
        components: AssetRegistryPathComponents,
        digest_on: "Callable[[str], GenomeDigest | None]",
        local_digest: GenomeDigest | None,
        server_urls: list[str] | None,
    ) -> str:
        """
        Remote seek, given how to find the genome's digest on each server.

        Args:
            components: The asset group, asset and seek key to find.
            digest_on: The genome's digest on one server (by URL), or None if
                that server cannot identify the genome.
            local_digest: The genome's digest, if the caller knows it here. A
                local non-path seek key value answers without a server.
            server_urls: Servers to query. Defaults to the subscriptions.
        """
        if local_digest is not None:
            local_value = self._local_non_path_value(local_digest, components)
            if local_value is not None:
                return local_value

        server_urls = server_urls or self.subscriptions()
        for server_url in server_urls:
            genome_digest = digest_on(server_url)
            if genome_digest is None:
                continue
            client = self.client(server_url)

            try:
                asset_groups = client.get_asset_groups(
                    params={
                        "genome_digest": genome_digest,
                        "asset_group_name": components.asset_group,
                    }
                )
            except Exception:
                continue
            if not asset_groups:
                continue

            # Fetch the group's assets from the server and choose one from the
            # response, rather than consulting the local DB (which need not have
            # this genome).
            asset_group_id = asset_groups[0]["id"]
            try:
                if components.asset is not None:
                    assets = client.get_assets(
                        params={"asset_group_id": asset_group_id, "name": components.asset}
                    )
                else:
                    assets = client.get_assets(params={"asset_group_id": asset_group_id})
            except Exception:
                continue
            if not assets:
                continue

            if components.asset is not None:
                asset_metadata = assets[0]
            else:
                # Default asset: the one flagged is_default, else the sole asset.
                default_assets = [a for a in assets if a.get("is_default") is True]
                if len(default_assets) == 1:
                    asset_metadata = default_assets[0]
                elif len(assets) == 1:
                    asset_metadata = assets[0]
                else:
                    continue  # ambiguous: no default and more than one asset

            serving_modes = asset_metadata.get("serving_modes", ["archive"])
            if "file" not in serving_modes:
                raise ValueError(
                    f"Asset '{components}' does not support file-level access "
                    f"(serving_modes={serving_modes}). Pull it locally with 'refgenie pull'."
                )

            # Choose the seek key from the server payload, not the local DB.
            seek_keys = asset_metadata.get("seek_keys") or []
            if components.seek_key is not None:
                chosen = next(
                    (sk for sk in seek_keys if sk.get("name") == components.seek_key),
                    None,
                )
            else:
                chosen = default_seek_key_from_payload(seek_keys, components.asset_group)
            if chosen is None:
                continue

            # Non-path seek keys return their value directly, no URL needed.
            if not payload_is_path_type(chosen.get("type")):
                return chosen["value"]

            return f"{server_url}/v4/assets/{asset_metadata['digest']}/files/{chosen['value']}"

        raise ValueError("Failed to resolve remote seek path from any subscribed server.")
