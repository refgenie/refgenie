"""
AssetAliasTreeMixin - the name-addressed alias/build trees of the AssetManager.

Owns rendering and purging of the directories an asset name owns per genome
alias: its alias-tree view (``alias/<alias>/<group>/<name>/``) and its build
bookkeeping (``builds/<alias>/<group>/<name>/``). Split out on the call graph:
both ``content.py``'s commit path and ``manager.py``'s remove/rename call in
here.

This is a mixin, not a collaborator object, so every method keeps its name on
``AssetManager``; it is not usable standalone.
"""

import shutil
from pathlib import Path
from typing import TYPE_CHECKING

from refgenie.exceptions import MissingAssetError
from refgenie.logger import logger

if TYPE_CHECKING:
    from refgenie.managers.alias import AliasManager


class AssetAliasTreeMixin:
    """
    Alias-tree renderer/purger mixed into ``AssetManager``.

    Relies on state and methods provided by ``AssetManager`` (attributes below,
    plus catalog lookups ``get`` and ``list_assets``), resolved through the MRO
    at runtime.
    """

    # Provided by AssetManager.__init__.
    _genome_folder: Path
    _alias_folder: Path
    _alias_manager: "AliasManager"

    def _render_asset_alias_symlinks(
        self,
        genome_digest: str,
        asset_group_name: str,
        asset_name: str,
        alias_whitelist: list[str] | None = None,
    ) -> None:
        """
        Render the alias-tree symlinks for one asset name across genome aliases.

        The alias tree is name-addressed: ``alias/<alias>/<group>/<name>/`` is a
        view of the content directory with filenames rewritten per alias. This
        renders that view for one (asset name) x every genome alias.
        """
        from refgenie.utils.symlinks import create_alias_symlinks, get_symlink_paths

        try:
            asset = self.get(
                genome_digest=genome_digest,
                asset_group_name=asset_group_name,
                asset_name=asset_name,
            )
        except MissingAssetError:
            return
        if asset.path is None:
            return
        src_path = self._genome_folder / asset.path
        aliases = [
            a
            for a in self._alias_manager.get_for_genome(genome_digest)
            if (alias_whitelist is None or a in alias_whitelist)
        ]
        if not aliases:
            return
        target_paths = get_symlink_paths(
            alias_folder=self._alias_folder,
            aliases=aliases,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
        )
        create_alias_symlinks(
            src_path=src_path,
            target_paths_mapping=target_paths,
            genome_digest=genome_digest,
        )

    def render_alias_tree(
        self,
        genome_digest: str,
        alias_whitelist: list[str] | None = None,
    ) -> None:
        """
        Render the full name-addressed alias tree for a genome.

        Used when a new alias is attached to an existing genome: the content
        lives under digest-named directories, so the alias view must be rebuilt
        from the ``assetname`` rows rather than by walking the data directory
        (which would produce digest-named alias directories).
        """
        for asset in self.list_assets(genome_digests=[genome_digest]):
            if asset.path is None:
                continue
            group_name = asset.asset_group.name
            for asset_name in asset.names:
                self._render_asset_alias_symlinks(
                    genome_digest, group_name, asset_name.name, alias_whitelist=alias_whitelist
                )

    def _owned_trees(
        self,
        genome_digest: str,
        asset_group_name: str,
        asset_name: str | None = None,
    ) -> list[tuple[str, Path]]:
        """
        The name-addressed directories owned by one asset name, or by a group.

        An asset name owns two directories per genome alias: its alias-tree view
        (``alias/<alias>/<group>/<name>/``) and its build bookkeeping
        (``builds/<alias>/<group>/<name>/``). Passing ``asset_name=None``
        addresses the group-level parents of both instead.

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

    def _purge_name_trees(
        self,
        genome_digest: str,
        asset_group_name: str,
        names: list[str],
    ) -> None:
        """
        Delete the alias-tree and build-tree directories for each of ``names``.

        Removal destroys the content and every name that pointed at it, and each
        of those names has its own pair of directories.

        The build directory holds the completion flag, which is not a historical
        record: it asserts that the asset currently exists. Its lifetime is
        deliberately tied to the asset's. Leave it behind and snakemake sees its
        declared output present, skips the rebuild, and the removed asset
        silently never comes back.
        """
        for name in names:
            for alias_name, path in self._owned_trees(genome_digest, asset_group_name, name):
                logger.info(f"Removing files for alias '{alias_name}': {path}")
                shutil.rmtree(path, ignore_errors=True)
