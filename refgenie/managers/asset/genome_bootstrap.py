"""
Making a genome exist locally before its assets are pulled into it.

A pull targets ``genome/asset_group:asset``, but the genome row and its alias
may not be here yet -- the server is the authority for both. These two methods
resolve the alias against the server, create what is missing, and report what
they created so `PullTransaction` can roll it back if the pull then fails.

Mixed into `AssetPuller`, so every method keeps its name and no call site
changes.
"""

from dataclasses import dataclass

from refgenie.exceptions import MissingGenomeError
from refgenie.logger import logger
from refgenie.managers.sources import make_source
from refgenie.managers.sources.api_ids import (
    API_ID_ALIAS_DIGEST,
    API_ID_GENOME_ATTRS,
)


@dataclass
class GenomeCreationResult:
    """Result of _ensure_genome_exists() with creation tracking."""

    success: bool
    genome_digest: str | None = None
    alias_name: str | None = None
    created_genome: bool = False
    created_alias: bool = False


class GenomeBootstrapMixin:
    """Genome/alias creation from a remote server. Mixed into `AssetPuller`."""

    def _ensure_genome_exists(
        self,
        alias_name: str,
        genome_digest: str | None = None,
        genome_description: str | None = None,
        server_urls: list[str] | None = None,
    ) -> GenomeCreationResult:
        """
        Ensure a genome exists for the given alias, creating it if necessary.

        This is used during pull operations. Unlike Refgenie.set_genome_alias(),
        this does NOT create symlinks - those are handled by _create_symlinks_for_alias()
        after asset extraction.

        Uses RemoteGenomeSource + GenomeManager.initialize_genome() so that
        remote_url is set on the Genome record (enabling lazy sequence fetching).

        Args:
            alias_name: The name of the alias.
            genome_digest: The digest of the genome. Optional.
            genome_description: The description of the genome.
            server_urls: The URLs of the server.

        Returns:
            GenomeCreationResult: Result with success flag and creation tracking info.
        """
        servers = server_urls or self._sources.get_subscriptions()

        if genome_digest is None:
            # Need to resolve alias via remote servers
            resolved_digest = None
            genome_description = genome_description or ""

            for server_url in servers:
                try:
                    source = make_source(server_url)
                except Exception as exc:
                    logger.debug(f"Failed to create source for {server_url}: {exc}")
                    continue
                resolved_digest = source.resolve_alias(alias_name)
                if resolved_digest is not None:
                    genome_digest = resolved_digest
                    logger.info(f"Resolved alias '{alias_name}' to digest: {genome_digest}")
                    break
                else:
                    logger.debug(f"Alias '{alias_name}' not found on {server_url}")

            if genome_digest is None:
                # Fall back to server client for alias resolution (v4 API)
                for server_url in servers:
                    client = self._sources.get_server_client(server_url)
                    try:
                        genome = client.get(
                            operation_id=API_ID_ALIAS_DIGEST,
                            url_format_params={"name": alias_name},
                        )
                        genome_digest = genome["digest"]
                        genome_attrs = client.get(
                            operation_id=API_ID_GENOME_ATTRS,
                            url_format_params={"digest": genome_digest},
                        )
                        if isinstance(genome_attrs, dict) and "description" in genome_attrs:
                            genome_description = genome_attrs["description"]
                    except ValueError as exc:
                        logger.debug(f"Server {server_url} does not support alias endpoint: {exc}")
                        continue
                    except Exception as exc:
                        logger.debug(
                            f"Alias resolution via server client failed for {server_url}: {exc}"
                        )
                        continue
                    if genome_digest:
                        logger.info(
                            f"Resolved alias '{alias_name}' via server client: {genome_digest}"
                        )
                        break

            if genome_digest is None:
                return GenomeCreationResult(success=False)

            logger.info(f"Determined digest for {alias_name}: {genome_digest}")

        # Record what already exists BEFORE any writes: the result must report
        # only rows this call actually created, so a pull rollback never
        # deletes a genome or alias that pre-existed the pull (the digest may
        # already be registered locally under a different alias).
        genome_existed = self._genome_manager.exists(genome_digest)
        alias_existed = self._alias_manager.exists(alias_name)

        try:
            self._genome_manager.get(genome_digest)
        except MissingGenomeError:
            if genome_description is None:
                genome_description = "No description provided"

            # Try to initialize via RemoteGenomeSource for remote_url support
            source = None
            for server_url in servers:
                try:
                    source = make_source(server_url)
                    break
                except Exception:
                    continue

            if source is not None:
                try:
                    _digest, created = self._genome_manager.initialize_genome(
                        source=source,
                        digest=genome_digest,
                        description=genome_description,
                        alias_names=[alias_name],
                        use_existing=True,
                    )
                    return GenomeCreationResult(
                        success=True,
                        genome_digest=genome_digest,
                        alias_name=alias_name,
                        created_genome=created,
                        created_alias=not alias_existed,
                    )
                except Exception as exc:
                    logger.debug(f"Source-based genome init failed, falling back: {exc}")

            # Fallback: direct add without remote_url
            self._genome_manager.add(genome_digest, genome_description, [alias_name])
            return GenomeCreationResult(
                success=True,
                genome_digest=genome_digest,
                alias_name=alias_name,
                created_genome=not genome_existed,
                created_alias=not alias_existed,
            )
        else:
            # Genome exists, add alias to it
            self._alias_manager.add(alias_name, genome_digest)
            # Note: Symlinks are created later in _create_symlinks_for_alias() after asset extraction
            return GenomeCreationResult(
                success=True,
                genome_digest=genome_digest,
                alias_name=alias_name,
                created_genome=not genome_existed,
                created_alias=not alias_existed,
            )

    def init_genome_from_remote(
        self,
        alias_name: str,
        genome_digest: str | None = None,
        genome_description: str | None = None,
        server_urls: list[str] | None = None,
    ) -> bool:
        """
        Register a genome in the local store from remote metadata.

        Creates genome entry and alias without downloading any asset files.
        Uses RemoteGenomeSource for remote_url support.

        Args:
            alias_name: The alias name to look up on remote servers.
            genome_digest: The digest of the genome, if already known.
            genome_description: The description of the genome.
            server_urls: Optional list of server URLs to query.

        Returns:
            True if genome was registered, False if not found or already exists.
        """
        if self._alias_manager.exists(alias_name):
            logger.info(f"Alias '{alias_name}' already exists locally")
            return True

        result = self._ensure_genome_exists(
            alias_name=alias_name,
            genome_digest=genome_digest,
            genome_description=genome_description,
            server_urls=server_urls,
        )
        if result.success:
            logger.info(f"Registered genome for alias '{alias_name}': {result.genome_digest}")
            return True
        else:
            logger.error(f"Could not find alias '{alias_name}' on remote servers")
            return False
