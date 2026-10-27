"""
Tests for the per-alias ``alias/`` and ``builds/`` trees (refgenie.managers.asset.alias_tree).

``AliasTree`` is reached as ``rgc.asset.tree``. One test builds it without an
AssetManager, rendering into a folder of its own, which proves its constructor
names everything it needs. The ``Refgenie``-level behavior around the tree
(``add`` rendering it, ``set_genome_alias`` refreshing it) stays in
test_genome.py; build-directory lookup stays in test_build_dir_location.py.
"""

from pathlib import Path

import pytest

from refgenie import Refgenie
from refgenie.managers.asset.alias_tree import AliasTree
from refgenie.managers.asset.group import AssetGroupManager
from tests.helpers import ASSET, GENOME, GROUP, add_asset_from_files, genome_digest

# Component tier: real genome folders and asset files.
pytestmark = pytest.mark.component

SECOND_ALIAS = "rCRSd_alt"
FASTA_FILES = ["{genome}.fa", "{genome}.fa.fai", "{genome}.chrom.sizes"]


def _alias_dir(r: Refgenie, alias: str, group: str, asset: str) -> Path:
    return r.alias_folder / alias / group / asset


def _add_fasta(r: Refgenie, asset_name: str = ASSET) -> None:
    add_asset_from_files(
        r,
        asset_class_name="fasta",
        asset_group_name=GROUP,
        asset_name=asset_name,
        rel_dir=f"data/fasta_payload_{asset_name}",
        files=FASTA_FILES,
    )


class TestAliasTree:
    """Rendering the name-addressed alias view of digest-addressed content."""

    def test_works_without_an_asset_manager(self, refgenie_fs, tmp_path):
        """Built from its constructor arguments alone, it renders and enumerates."""
        r = refgenie_fs
        _add_fasta(r)
        digest = genome_digest(r)
        alias_folder = tmp_path / "separate_alias_folder"
        tree = AliasTree(
            r.database_engine,
            r.genome_folder,
            alias_folder,
            r.alias,
            AssetGroupManager(r.database_engine, r.genome),
        )

        tree.render(digest, GROUP, ASSET)

        view = alias_folder / GENOME / GROUP / ASSET
        assert (view / f"{GENOME}.fa").is_symlink()
        assert (GENOME, view) in tree.owned_trees(digest, GROUP, ASSET)

    def test_each_alias_tree_is_named_for_its_own_alias(self, refgenie_fs):
        """Regression test: a second alias must not inherit the first alias's names."""
        r = refgenie_fs
        digest = genome_digest(r)
        r.alias.add(name=SECOND_ALIAS, genome_digest=digest)
        _add_fasta(r)
        r.asset.tree.render(digest, GROUP, ASSET)

        for alias in (GENOME, SECOND_ALIAS):
            alias_dir = _alias_dir(r, alias, GROUP, ASSET)
            assert (alias_dir / f"{alias}.fa").is_symlink(), f"missing payload for {alias}"
            assert not (alias_dir / f"{digest}.fa").exists(), f"digest survived in {alias} tree"

        # The decisive assertion: neither tree is named for the other's alias.
        assert not (_alias_dir(r, SECOND_ALIAS, GROUP, ASSET) / f"{GENOME}.fa").exists()

    def test_render_without_an_asset_name_uses_the_group_default(self, refgenie_fs):
        """Rendering a group with no asset name resolves the group's default
        asset (by name, not digest) and renders its tree."""
        r = refgenie_fs
        digest = genome_digest(r)
        _add_fasta(r)
        assert r.asset.group.get_default(GROUP, genome_digest=digest) == ASSET

        r.asset.tree.render(digest, GROUP)

        assert _alias_dir(r, GENOME, GROUP, ASSET).is_dir()

    def test_render_of_a_missing_asset_is_a_no_op(self, refgenie_fs):
        """An unknown asset name renders nothing and raises nothing."""
        r = refgenie_fs
        digest = genome_digest(r)
        _add_fasta(r)
        r.asset.tree.render(digest, GROUP, "no_such_asset")
        assert not _alias_dir(r, GENOME, GROUP, "no_such_asset").exists()

    def test_render_on_an_alias_less_genome_warns_and_touches_nothing(self, refgenie_fs, caplog):
        """A genome with no aliases has nowhere to render the tree: warn, don't
        raise, and don't read the catalog for the named asset."""
        r = refgenie_fs
        digest = genome_digest(r)
        _add_fasta(r)
        r.alias.remove(GENOME)

        with caplog.at_level("WARNING"):
            r.asset.tree.render(digest, GROUP, ASSET)

        assert any("has no aliases; skipping alias tree" in rec.message for rec in caplog.records)
        assert not (r.alias_folder / GENOME / GROUP / ASSET).exists()

    def test_render_genome_honors_the_alias_whitelist(self, refgenie_fs):
        """render_genome renders every name of every asset, for only the listed aliases."""
        r = refgenie_fs
        digest = genome_digest(r)
        _add_fasta(r)
        r.alias.add(name=SECOND_ALIAS, genome_digest=digest)
        assert not _alias_dir(r, SECOND_ALIAS, GROUP, ASSET).exists()

        r.asset.tree.render_genome(digest, aliases=[SECOND_ALIAS])

        assert (_alias_dir(r, SECOND_ALIAS, GROUP, ASSET) / f"{SECOND_ALIAS}.fa").is_symlink()
