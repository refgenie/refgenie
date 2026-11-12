"""
AliasTree - the name-addressed alias and build trees on disk.

Reached as ``rgc.asset.tree``. Owns rendering and purging of the directories an
asset name owns per genome alias: its alias-tree view
(``alias/<alias>/<group>/<name>/``) and its build bookkeeping
(``builds/<alias>/<group>/<name>/``). It writes symlinks and directories, not
catalog rows, and reads the catalog directly rather than through
``AssetManager``.
"""

import shutil
from pathlib import Path
from typing import TYPE_CHECKING

from refgenie.db.engine import Engine
from refgenie.db.queries import asset as queries
from refgenie.logger import logger
from refgenie.managers.base import ResourceManager
from refgenie.models import GenomeAlias, GenomeDigest

if TYPE_CHECKING:
    from refgenie.managers.alias import AliasManager
    from refgenie.managers.asset.group import AssetGroupManager


class AliasTree(ResourceManager):
    """Renderer and purger of the per-alias ``alias/`` and ``builds/`` trees."""

    def __init__(
        self,
        database_engine: Engine,
        genome_folder: Path,
        alias_folder: Path,
        alias_manager: "AliasManager",
        groups: "AssetGroupManager",
    ):
        """
        Initialize the AliasTree.

        Args:
            database_engine: The database engine.
            genome_folder: Path to genome data folder (holds ``data/`` and ``builds/``).
            alias_folder: Path to alias folder.
            alias_manager: The AliasManager, for the aliases of a genome.
            groups: The AssetGroupManager, for a group's default asset.
        """
        super().__init__(database_engine)
        self._genome_folder = genome_folder
        self._alias_folder = alias_folder
        self._alias_manager = alias_manager
        self._groups = groups

    def _aliases_for(
        self, genome_digest: GenomeDigest, aliases: list[GenomeAlias] | None = None
    ) -> list[str]:
        """
        Get the aliases for a genome, optionally restricted to a whitelist.

        Args:
            genome_digest: The digest of the genome.
            aliases: The list of aliases to include.

        Returns:
            list[str]: The matching alias names.
        """
        return [
            alias
            for alias in self._alias_manager.get_for_genome(genome_digest)
            if (aliases is None or alias in aliases)
        ]

    def _link(
        self,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        asset_name: str,
        asset_path: str,
        aliases: list[GenomeAlias] | None,
    ) -> None:
        """Render one asset name's view of ``asset_path`` under each selected alias."""
        from refgenie.utils.symlinks import create_alias_symlinks, get_symlink_paths

        alias_names = self._aliases_for(genome_digest, aliases)
        if not alias_names:
            return
        target_paths = get_symlink_paths(
            alias_folder=self._alias_folder,
            aliases=alias_names,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
        )
        create_alias_symlinks(
            src_path=self._genome_folder / asset_path,
            target_paths_mapping=target_paths,
            genome_digest=genome_digest,
        )

    def render(
        self,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        asset_name: str | None = None,
    ) -> None:
        """
        Render the alias-tree symlinks for one asset name across genome aliases.

        ``alias/<alias>/<group>/<name>/`` is a view of the content directory
        with filenames rewritten per alias. With no ``asset_name`` the group's
        default is rendered. Logs a warning if there is no default or the
        genome has no aliases; returns quietly if the asset is missing.

        Raises:
            MissingAssetGroupError: If ``asset_name`` is omitted and the group
                does not exist.
        """
        asset_name = asset_name or self._groups.get_default(
            asset_group_name, genome_digest=genome_digest
        )
        if asset_name is None:
            logger.warning(
                f"Not rendering the alias tree for '{genome_digest}/{asset_group_name}': "
                f"no asset name was given and the group has no default asset."
            )
            return
        if not self._alias_manager.get_for_genome(genome_digest):
            logger.warning(
                f"Genome '{genome_digest}' has no aliases; skipping alias tree for "
                f"'{asset_group_name}:{asset_name}'"
            )
            return
        with self._database_session as session:
            asset = queries.find_by_name(session, genome_digest, asset_group_name, asset_name)
        if asset is None:
            return
        self._link(genome_digest, asset_group_name, asset_name, asset.path, None)

    def render_genome(
        self, genome_digest: GenomeDigest, aliases: list[GenomeAlias] | None = None
    ) -> None:
        """
        Render the full name-addressed alias tree for a genome.

        Used when a new alias is attached to an existing genome. Built from the
        ``assetname`` rows, since the data directory is digest-named.
        """
        with self._database_session as session:
            assets = queries.list_all(session, [genome_digest])
        for asset in assets:
            group_name = asset.asset_group.name
            for asset_name in asset.names:
                self._link(genome_digest, group_name, asset_name.name, asset.path, aliases)

    def build_paths(
        self,
        genome_digest: GenomeDigest,
        asset_group_name: str | None = None,
        asset_name: str | None = None,
        aliases: list[GenomeAlias] | None = None,
    ) -> dict[str, Path]:
        """
        Get path to the build directory for the selected genome-group-asset.

        Mirrors :func:`refgenie.utils.symlinks.get_symlink_paths`, but rooted at
        the ``builds/`` tree.

        Args:
            genome_digest: The digest of the genome.
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            aliases: The list of aliases to include.

        Returns:
            dict[str, Path]: Mapping of alias names to their build directory paths.
        """
        from refgenie.utils.symlinks import get_build_paths

        return get_build_paths(
            genome_folder=self._genome_folder,
            aliases=self._aliases_for(genome_digest, aliases),
            asset_group_name=asset_group_name,
            asset_name=asset_name,
        )

    def find_build_dir(
        self,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        asset_name: str,
    ) -> Path | None:
        """
        Find the build directory for an asset, whichever alias it was built under.

        ``builds/`` is keyed by the alias used at build time, which may differ
        from the one a later command names, so every alias is searched.

        Args:
            genome_digest: The digest of the genome.
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.

        Returns:
            Path | None: The build directory, or None if there is none (e.g. a
            pulled asset).
        """
        for path in self.build_paths(
            genome_digest=genome_digest,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
        ).values():
            if path.is_dir():
                return path
        return None

    def owned_trees(
        self,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        asset_name: str | None = None,
    ) -> list[tuple[str, Path]]:
        """
        The name-addressed directories owned by one asset name, or by a group.

        Per genome alias, a name owns ``alias/<alias>/<group>/<name>/`` and
        ``builds/<alias>/<group>/<name>/``. ``asset_name=None`` gives the
        group-level parents instead.

        Returns:
            list[tuple[str, Path]]: (alias name, directory) pairs.
        """
        from refgenie.utils.symlinks import get_build_paths, get_symlink_paths

        aliases = self._alias_manager.get_for_genome(genome_digest)
        symlink_paths = get_symlink_paths(
            alias_folder=self._alias_folder,
            aliases=aliases,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
        )
        build_paths = get_build_paths(
            genome_folder=self._genome_folder,
            aliases=aliases,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
        )
        return list(symlink_paths.items()) + list(build_paths.items())

    def purge_names(
        self,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        names: list[str],
    ) -> None:
        """
        Delete the alias-tree and build-tree directories for each of ``names``.

        The build directory holds the completion flag, so it must go with the
        asset; otherwise snakemake would skip a later rebuild.
        """
        for name in names:
            for alias_name, path in self.owned_trees(genome_digest, asset_group_name, name):
                logger.info(f"Removing files for alias '{alias_name}': {path}")
                shutil.rmtree(path, ignore_errors=True)
