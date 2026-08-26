"""
Reading a remote server's catalog: what genomes and assets it has, and how big.

Split out of `AssetPuller` because none of it touches the filesystem or a
database transaction -- it queries subscribed servers and shapes the answers
into dicts. Mixed back into `AssetPuller`, so every method keeps its name and
no call site changes.
"""

from collections import defaultdict

from refgenie.logger import logger


class AssetRemoteCatalogMixin:
    """Server-side catalog queries. Mixed into `AssetPuller`."""

    def list_remote(
        self,
        genome_digests: list[str] | None = None,
        include_seek_keys: bool = False,
        server_urls: list[str] | None = None,
    ) -> tuple[dict[str, dict[str, list[str]]], dict[str, dict[str, str]]]:
        """
        List all assets available on the server subscriptions.

        Args:
            genome_digests: Optional list of genome digests to filter by.
            include_seek_keys: Whether to include seek keys in the results.
            server_urls: Optional list of server URLs to query. If not provided,
                         uses all subscribed servers.

        Returns:
            tuple[dict[str, dict[str, list[str]]], dict[str, dict[str, str]]]: A tuple containing:
                - Dictionary with server URLs as keys and asset data as values
                - Dictionary with server URLs as keys and alias data as values
                Asset data follows the same format as the list() method:
                genome digests as keys and lists of asset strings as values.
        """
        server_urls = server_urls or self._sources._config.get_server_subscriptions()
        if not server_urls:
            logger.warning("No server subscriptions found")
            return {}, {}

        result = {}
        aliases_result = {}

        for server_url in server_urls:
            logger.debug(f"Querying remote assets from: {server_url}")
            client = self.get_client(server_url)

            # Query ALL assets from the server using pagination
            assets_response = client.get_all_assets()
            if not assets_response:
                logger.warning(f"No assets found on server: {server_url}")
                result[server_url] = {}
                aliases_result[server_url] = {}
                continue

            # Query ALL asset groups from the server using pagination
            asset_groups_response = client.get_all_asset_groups()

            # Query ALL genomes from the server using pagination
            genomes_response = client.get_all_genomes()

            # Query ALL aliases from the server using pagination
            try:
                aliases_response = client.get_all_aliases()
            except Exception as e:
                logger.warning(f"Failed to fetch aliases from {server_url}: {e}")
                aliases_response = []

            # Build lookup dictionaries for in-memory matching
            asset_groups = {}
            for asset_group in asset_groups_response or []:
                asset_groups[asset_group["id"]] = asset_group

            genomes = {}
            for genome in genomes_response or []:
                genomes[genome["digest"]] = genome

            # Group assets by asset_group_id and filter by genome_digests if provided
            assets_by_group = {}
            for asset in assets_response:
                asset_group_id = asset.get("asset_group_id")
                if not asset_group_id:
                    continue

                # Get the asset group to check genome_digest
                asset_group = asset_groups.get(asset_group_id)
                if not asset_group:
                    continue

                asset_genome_digest = asset_group.get("genome_digest")

                # Filter by genome_digests if provided
                if genome_digests and asset_genome_digest not in genome_digests:
                    continue

                if asset_group_id not in assets_by_group:
                    assets_by_group[asset_group_id] = []
                assets_by_group[asset_group_id].append(asset)

            # Build the result structure similar to local list() method
            server_assets = {}
            for asset_group_id, group_assets in assets_by_group.items():
                asset_group = asset_groups.get(asset_group_id, {})
                genome_digest = asset_group.get("genome_digest")
                asset_group_name = asset_group.get("name", "unknown")

                if genome_digest:
                    if genome_digest not in server_assets:
                        server_assets[genome_digest] = []

                    for asset in group_assets:
                        asset_name = asset.get("name", "unknown")
                        if include_seek_keys:
                            # Note: seek keys are not currently available from remote APIs
                            # so we default to the same format as local list() method
                            server_assets[genome_digest].append(f"{asset_group_name}:{asset_name}")
                        else:
                            server_assets[genome_digest].append(f"{asset_group_name}:{asset_name}")

            # Build aliases lookup for this server
            server_aliases = {}
            genome_digests_with_assets = list(server_assets.keys())
            # Group aliases by genome digest
            aliases_by_genome = defaultdict(list)
            for alias_item in aliases_response:
                genome_digest = alias_item.get("genome_digest")
                alias_name = alias_item.get("name")
                if genome_digest in genome_digests_with_assets and alias_name:
                    aliases_by_genome[genome_digest].append(alias_name)

            # Convert to comma-separated strings
            for genome_digest, alias_names in aliases_by_genome.items():
                server_aliases[genome_digest] = ", ".join(alias_names)

            result[server_url] = server_assets
            aliases_result[server_url] = server_aliases
            logger.debug(f"Found {len(server_assets)} genomes with assets on {server_url}")

        return result, aliases_result

    def list_remote_assets_for_genome(
        self,
        genome_digest: str,
        server_urls: list[str] | None = None,
    ) -> list[dict]:
        """
        List all assets available for a genome on remote servers.

        Args:
            genome_digest: The digest of the genome to query.
            server_urls: Optional list of server URLs to query. If not provided,
                         uses all subscribed servers.

        Returns:
            List of dicts with keys: server_url, genome_digest, asset_group_name,
            asset_name, archive_size, asset_digest, archive_digest.
        """
        server_urls = server_urls or self._sources.get_subscriptions()
        if not server_urls:
            logger.warning("No server subscriptions found")
            return []

        results = []
        for server_url in server_urls:
            logger.debug(f"Querying remote assets for genome {genome_digest} from: {server_url}")
            client = self.get_client(server_url)

            # Get asset groups for this genome
            try:
                asset_groups = client.get_all_asset_groups(params={"genome_digest": genome_digest})
            except Exception as e:
                logger.warning(f"Failed to fetch asset groups from {server_url}: {e}")
                continue

            if not asset_groups:
                continue

            for asset_group in asset_groups:
                asset_group_id = asset_group.get("id")
                asset_group_name = asset_group.get("name", "unknown")

                # Get assets for this asset group
                try:
                    assets = client.get_all_assets(params={"asset_group_id": asset_group_id})
                except Exception as e:
                    logger.warning(f"Failed to fetch assets from {server_url}: {e}")
                    continue

                for asset in assets or []:
                    asset_name = asset.get("name", "unknown")
                    asset_digest = asset.get("digest")

                    # Get archive size for this asset
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
                            "asset_name": asset_name,
                            "asset_digest": asset_digest,
                            "archive_digest": archive_digest,
                            "archive_size": archive_size,
                        }
                    )

        return results

    def list_remote_genomes(
        self,
        server_urls: list[str] | None = None,
    ) -> list[dict]:
        """
        List all genomes available on remote servers.

        Args:
            server_urls: Optional list of server URLs to query. If not provided,
                         uses all subscribed servers.

        Returns:
            List of dicts with keys: server_url, genome_digest, aliases, description.
        """
        server_urls = server_urls or self._sources.get_subscriptions()
        if not server_urls:
            logger.warning("No server subscriptions found")
            return []

        results = []
        seen_digests = set()

        for server_url in server_urls:
            logger.debug(f"Querying remote genomes from: {server_url}")
            client = self.get_client(server_url)

            try:
                genomes = client.get_all_genomes()
            except Exception as e:
                logger.warning(f"Failed to fetch genomes from {server_url}: {e}")
                continue

            try:
                aliases = client.get_all_aliases()
            except Exception as e:
                logger.warning(f"Failed to fetch aliases from {server_url}: {e}")
                aliases = []

            # Build alias lookup by genome digest
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

    def estimate_pull_size(
        self,
        asset_list: list[dict],
    ) -> tuple[int, int]:
        """
        Calculate total size and count for a list of remote assets.

        Args:
            asset_list: List of asset dicts from list_remote_assets_for_genome().
                        Each dict should have an 'archive_size' key.

        Returns:
            Tuple of (total_bytes, asset_count). If archive_size is unknown for
            some assets, they are counted but their size is not included in the total.
        """
        total_bytes = 0
        asset_count = len(asset_list)

        for asset in asset_list:
            size = asset.get("archive_size")
            if size is not None:
                total_bytes += size

        return total_bytes, asset_count
