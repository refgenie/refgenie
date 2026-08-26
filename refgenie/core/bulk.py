"""Bulk transfer operations: pull-all, per-asset bulk pulls, genome registration, mirror."""

from refgenie.logger import logger
from refgenie.utils.prompt import Confirmer


class BulkTransferMixin:
    """
    Bulk pull/mirror operations over many genomes and assets.

    Mixed into :class:`refgenie.core.facade.Refgenie`; relies on the facade's
    ``alias`` and ``asset`` manager properties.
    """

    def pull_all(
        self,
        genome_names: list[str],
        force: bool | None = None,
        force_large: bool | None = None,
        size_cutoff: int | float | None = None,
        confirm: Confirmer | None = None,
    ) -> list:
        """
        Pull all assets for the specified genome(s).

        Shows size estimate and confirmation prompt (unless force=True).

        Args:
            genome_names: List of genome aliases or digests.
            force: Skip confirmation prompts if True.
            force_large: How to handle large archives.
            size_cutoff: Maximum archive file size to download without prompt.

        Returns:
            List of successfully pulled Assets.
        """
        from refgenie.utils.io import confirm_bulk_pull

        genome_digests = []
        for name in genome_names:
            try:
                digest = self.alias.resolve(name)
            except Exception:
                # Try resolving from remote
                result = self.asset.init_genome_from_remote(name)
                if result:
                    digest = self.alias.resolve(name)
                else:
                    logger.error(f"Could not resolve genome: {name}")
                    continue
            genome_digests.append(digest)

        if not genome_digests:
            logger.error("No valid genomes to pull")
            return []

        all_assets = []
        for digest in genome_digests:
            assets = self.asset.list_remote_assets_for_genome(digest)
            all_assets.extend(assets)

        if not all_assets:
            logger.warning("No remote assets found for specified genomes")
            return []

        total_bytes, asset_count = self.asset.estimate_pull_size(all_assets)
        if not confirm_bulk_pull(
            asset_count=asset_count,
            genome_count=len(genome_digests),
            total_bytes=total_bytes,
            force=force or False,
            confirm=confirm,
        ):
            logger.info("Bulk pull cancelled")
            return []

        return self.asset.pull_multiple(
            asset_list=all_assets,
            force=force,
            force_large=force_large,
            size_cutoff=size_cutoff,
            confirm=confirm,
        )

    def pull_asset_for_genomes(
        self,
        asset_name: str,
        genome_names: list[str] | None = None,
        all_genomes: bool = False,
        force: bool | None = None,
        force_large: bool | None = None,
        size_cutoff: int | float | None = None,
        confirm: Confirmer | None = None,
    ) -> list:
        """
        Pull a specific asset for multiple genomes.

        Args:
            asset_name: Name of the asset to pull (e.g., "fasta").
            genome_names: List of genome aliases or digests. Ignored if all_genomes=True.
            all_genomes: If True, queries servers for all available genomes.
            force: Skip confirmation prompts if True.
            force_large: How to handle large archives.
            size_cutoff: Maximum archive file size to download without prompt.

        Returns:
            List of successfully pulled Assets.
        """
        from refgenie.utils.io import confirm_bulk_pull

        if all_genomes:
            remote_genomes = self.asset.list_remote_genomes()
            genome_names = [
                g["aliases"][0] if g["aliases"] else g["genome_digest"] for g in remote_genomes
            ]

        if not genome_names:
            logger.error("No genomes specified")
            return []

        genome_digests = []
        for name in genome_names:
            try:
                digest = self.alias.resolve(name)
            except Exception:
                result = self.asset.init_genome_from_remote(name)
                if result:
                    digest = self.alias.resolve(name)
                else:
                    logger.error(f"Could not resolve genome: {name}")
                    continue
            genome_digests.append(digest)

        if not genome_digests:
            logger.error("No valid genomes to pull")
            return []

        matching_assets = []
        for digest in genome_digests:
            assets = self.asset.list_remote_assets_for_genome(digest)
            for asset in assets:
                if asset.get("asset_group_name") == asset_name:
                    matching_assets.append(asset)

        if not matching_assets:
            logger.warning(f"No remote assets named '{asset_name}' found for specified genomes")
            return []

        total_bytes, asset_count = self.asset.estimate_pull_size(matching_assets)
        if not confirm_bulk_pull(
            asset_count=asset_count,
            genome_count=len(genome_digests),
            total_bytes=total_bytes,
            force=force or False,
            confirm=confirm,
        ):
            logger.info("Bulk pull cancelled")
            return []

        return self.asset.pull_multiple(
            confirm=confirm,
            asset_list=matching_assets,
            force=force,
            force_large=force_large,
            size_cutoff=size_cutoff,
        )

    def init_genomes(
        self,
        genome_names: list[str] | None = None,
        all_genomes: bool = False,
    ) -> list[str]:
        """
        Register genome(s) in local store from remote metadata.

        No files downloaded. Returns list of registered genome digests.

        Args:
            genome_names: List of genome aliases to register.
            all_genomes: If True, register all genomes from remote servers.

        Returns:
            List of registered genome digests.
        """
        registered = []

        if all_genomes:
            remote_genomes = self.asset.list_remote_genomes()
            for genome_info in remote_genomes:
                aliases = genome_info.get("aliases", [])
                if aliases:
                    alias = aliases[0]
                    if self.asset.init_genome_from_remote(alias):
                        try:
                            digest = self.alias.resolve(alias)
                            registered.append(digest)
                        except Exception:
                            pass
        elif genome_names:
            for name in genome_names:
                if self.asset.init_genome_from_remote(name):
                    try:
                        digest = self.alias.resolve(name)
                        registered.append(digest)
                    except Exception:
                        pass

        logger.info(f"Registered {len(registered)} genome(s)")
        return registered

    def mirror(
        self,
        force: bool | None = None,
        force_large: bool | None = None,
        size_cutoff: int | float | None = None,
        confirm: Confirmer | None = None,
    ) -> list:
        """
        Mirror all assets for all genomes from subscribed servers.

        Always shows confirmation with total size estimate.

        Args:
            force: Skip confirmation prompts if True.
            force_large: How to handle large archives.
            size_cutoff: Maximum archive file size to download without prompt.

        Returns:
            List of successfully pulled Assets.
        """
        from refgenie.utils.io import confirm_bulk_pull

        remote_genomes = self.asset.list_remote_genomes()
        if not remote_genomes:
            logger.warning("No remote genomes found")
            return []

        # Initialize all genomes first (metadata only)
        # Pass the digest that list_remote_genomes already retrieved
        for g in remote_genomes:
            alias_name = g["aliases"][0] if g.get("aliases") else None
            self.asset.init_genome_from_remote(
                alias_name=alias_name or g["genome_digest"],
                genome_digest=g["genome_digest"],
                genome_description=g.get("description"),
            )

        all_assets = []
        for g in remote_genomes:
            assets = self.asset.list_remote_assets_for_genome(g["genome_digest"])
            all_assets.extend(assets)

        if not all_assets:
            logger.warning("No remote assets found to mirror")
            return []

        total_bytes, asset_count = self.asset.estimate_pull_size(all_assets)
        if not confirm_bulk_pull(
            asset_count=asset_count,
            genome_count=len(remote_genomes),
            total_bytes=total_bytes,
            force=force or False,
            confirm=confirm,
        ):
            logger.info("Mirror cancelled")
            return []

        return self.asset.pull_multiple(
            asset_list=all_assets,
            force=force,
            force_large=force_large,
            size_cutoff=size_cutoff,
            confirm=confirm,
        )
