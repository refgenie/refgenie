"""
Asset and genome removal (component tier).

Removal is one algorithm reached through two entry points -- ``remove()`` (the
way the CLI calls it) and ``remove_by_digest()`` -- and it must tear down
everything a build created:

* the catalog rows (the asset, every one of its names, and the asset group
  once it empties -- a surviving group row is a group with no assets that
  nothing ever reaps);
* the digest-addressed content directory;
* the per-name alias directories in BOTH trees (``alias/`` and ``builds/``) --
  group-level cleanup previously substituted the group's stale default asset
  name, resolved to an already-deleted directory, and orphaned the group dirs;
* the build completion flag -- it asserts the asset currently exists; if it
  survives, snakemake sees its declared output present and the removed asset
  silently never rebuilds.

Tests assert on combined database + disk snapshots of one standard built
world, rather than one property per world-build.
"""

import pytest
from sqlmodel import Session, select

from refgenie import Refgenie
from refgenie.db.tables import Asset
from tests.helpers import (
    ASSET,
    GENOME,
    GROUP,
    boom,
    build_dir,
    build_flag,
    content_dir,
    db_snapshot,
    disk_snapshot,
    make_built_refgenie,
    only_asset,
)

# Component tier: real genome folders and real SQLite together, ~0.2-0.5s each.
# Deselected from the bare `pytest` inner loop; run with `pytest -m component`.
pytestmark = pytest.mark.component

SECOND_NAME = "second"


# --- The two removal entry points ----------------------------------------------


def _remove_by_name(r: Refgenie) -> None:
    r.asset.remove(
        genome_name=GENOME,
        genome_digest=r.alias.resolve(GENOME),
        asset_group_name=GROUP,
        asset_name=ASSET,
    )


def _remove_by_digest(r: Refgenie) -> None:
    r.asset.remove_by_digest(digest=only_asset(r).digest)


both_removal_paths = pytest.mark.parametrize(
    "remove", [_remove_by_name, _remove_by_digest], ids=["remove", "remove_by_digest"]
)


class TestLastAssetRemoval:
    """Removing a group's last asset tears down rows, content, aliases and flags."""

    @both_removal_paths
    def test_removal_clears_database_and_disk(self, refgenie_built, remove):
        """
        The combined teardown contract, for both entry points:

        * the asset, its name, and the emptied group row are deleted -- an
          orphaned assetgroup row is a group with no assets that nothing reaps;
        * the content directory is gone -- deferring removal past COMMIT must
          not turn it into never removing it;
        * the group directories are gone from both the alias tree and the
          builds tree -- the stale-default-name regression;
        * the completion flag is destroyed, or snakemake never rebuilds.
        """
        r = refgenie_built
        flag = build_flag(r)
        content = content_dir(r)
        alias_group_dir = r.alias_folder / GENOME / GROUP
        build_group_dir = build_dir(r, asset_name=None)

        # Preconditions: the built world has all of these.
        assert flag.is_file()
        assert content.is_dir()
        assert alias_group_dir.is_dir()
        assert build_group_dir.is_dir()

        remove(r)

        # The catalog: asset, name and group rows are gone, and the genome row
        # cascades away with its last group (the sequences live in the refget
        # store, not this table).
        assert db_snapshot(r) == {
            "genomes": [],
            "groups": [],
            "assets": [],
            "asset_names": [],
        }

        # The disk: content, both group directories and the flag are gone, and
        # no flag debris hides anywhere else under the genome folder.
        for gone in (flag, content, alias_group_dir, build_group_dir):
            assert not gone.exists(), f"survived removal: {gone}"
        leftover_flags = [p for p in disk_snapshot(r) if p.endswith(".flag")]
        assert leftover_flags == [], f"flag debris survived removal: {leftover_flags}"

    def test_remove_and_remove_by_digest_agree(self, tmp_path, fixtures_path):
        """
        The two removal entry points are one algorithm; they must not diverge.

        Two identical installations, one asset each, removed through the two
        different paths. Full database and on-disk state must match afterwards
        -- this catches divergence in anything the explicit assertions above
        do not name.
        """
        by_name = make_built_refgenie(tmp_path / "by_name", fixtures_path)
        by_digest = make_built_refgenie(tmp_path / "by_digest", fixtures_path)

        assert db_snapshot(by_name) == db_snapshot(by_digest)
        assert disk_snapshot(by_name) == disk_snapshot(by_digest)

        _remove_by_name(by_name)
        _remove_by_digest(by_digest)

        assert db_snapshot(by_name) == db_snapshot(by_digest)
        assert disk_snapshot(by_name) == disk_snapshot(by_digest)

    def test_remove_keeps_the_group_when_asked(self, refgenie_built):
        """
        ``keep_asset_group=True`` leaves the (now empty) group row in place --
        but the flag and build bookkeeping still die with the asset.
        """
        r = refgenie_built
        genome_digest = r.alias.resolve(GENOME)
        flag = build_flag(r)
        assert flag.is_file()

        r.asset.remove(
            genome_name=GENOME,
            genome_digest=genome_digest,
            asset_group_name=GROUP,
            asset_name=ASSET,
            keep_asset_group=True,
        )

        assert db_snapshot(r) == {
            "genomes": [genome_digest],
            "groups": [(genome_digest, GROUP)],
            "assets": [],
            "asset_names": [],
        }
        assert not flag.exists(), f"build flag survived asset removal: {flag}"
        assert not build_dir(r).exists()


class TestMultiNameRemoval:
    """One content digest can carry several names; removal must clear them all."""

    def test_remove_clears_every_alias_directory(self, refgenie_built):
        """
        A successful remove tears down the alias directory for every name.

        A rebuild under a new name produces byte-identical content, so it
        reuses the one Asset row and adds a second AssetName. Removing via the
        canonical name must clear the alias directories of both names, and
        leave no asset or name rows behind.
        """
        r = refgenie_built
        r.build_asset(
            recipe_name="fasta",
            genome_name=GENOME,
            asset_group_name=GROUP,
            asset_name=SECOND_NAME,
        )
        genome_digest = r.alias.resolve(GENOME)
        canonical = r.asset.get(
            genome_digest=genome_digest, asset_group_name=GROUP, asset_name=ASSET
        ).name

        r.asset.remove(
            genome_digest=genome_digest, asset_group_name=GROUP, asset_name=canonical
        )

        assert not (r.alias_folder / GENOME / GROUP / ASSET).exists()
        assert not (r.alias_folder / GENOME / GROUP / SECOND_NAME).exists()
        db = db_snapshot(r)
        assert db["assets"] == []
        assert db["asset_names"] == []


class TestDeleteOrdering:
    """The flush must not destroy content before COMMIT."""

    def test_rolled_back_delete_leaves_content_intact(self, refgenie_built):
        """
        ``Asset.before_delete`` used to ``shutil.rmtree`` the content directory
        during the flush -- before ``COMMIT``. ``session.rollback()`` restored
        the row and could not restore the bytes, leaving a catalog describing
        data that no longer existed.

        (The committed half of the contract -- a committed delete DOES remove
        the content -- is asserted by TestLastAssetRemoval above.)
        """
        r = refgenie_built
        content = content_dir(r)
        digest = only_asset(r).digest
        assert content.is_dir()

        with Session(r.database_engine) as session:
            row = session.exec(select(Asset).where(Asset.digest == digest)).unique().one()
            session.delete(row)
            session.flush()
            assert content.is_dir(), "the flush destroyed the content before COMMIT"
            session.rollback()

        assert content.is_dir(), "a rolled-back delete destroyed the content"
        assert r.asset.exists(
            genome_digest=r.alias.resolve(GENOME), asset_group_name=GROUP, asset_name=ASSET
        ), "the row did not come back"


class TestGenomeRemoval:
    """The catalog decides first; alias/file cleanup follows and is re-runnable."""

    def test_remove_commits_before_touching_the_aliases(self, refgenie_built, monkeypatch):
        """
        Aliases used to be removed before the genome row. An interruption between
        the two left every row intact but the genome unaddressable by name -- a live
        genome nobody could reach. The inverse leaves store aliases pointing at a
        genome that is gone, which no read path mistakes for a working genome.
        """
        r = refgenie_built
        digest = r.alias.resolve(GENOME)
        monkeypatch.setattr(r.alias, "remove", boom)

        with pytest.raises(RuntimeError):
            r.genome.remove(digest)

        assert not r.genome.exists(digest), (
            "the alias removal ran before the catalog committed its decision"
        )

    def test_remove_cleans_up_after_a_committed_decision(self, refgenie_built):
        """The reordering must not turn 'later' into 'never'."""
        r = refgenie_built
        digest = r.alias.resolve(GENOME)
        alias_dir = r.genome_folder / "alias" / GENOME
        assert alias_dir.exists()

        r.genome.remove(digest)

        assert not r.genome.exists(digest)
        assert db_snapshot(r)["genomes"] == []
        assert not alias_dir.exists(), "the alias tree outlived the genome"

    def test_removing_an_alias_twice_still_clears_its_trees(self, refgenie_built):
        """
        Killed between dropping the store alias and removing its files, a re-run
        used to raise ``MissingAliasError`` and leave the trees forever. The trees
        are derived from the name alone, so cleaning them is safe either way.
        """
        from refgenie.exceptions import MissingAliasError

        r = refgenie_built
        alias_dir = r.genome_folder / "alias" / GENOME
        build_tree = r.genome_folder / "builds" / GENOME
        assert alias_dir.exists() and build_tree.exists()

        # Drop the alias from the store only -- exactly the state a kill in the gap
        # leaves behind.
        r.refget_store.remove_collection_alias("refgenie", GENOME)
        r.alias.invalidate()
        assert alias_dir.exists(), "precondition: the trees are still there"

        with pytest.raises(MissingAliasError):
            r.alias.remove(GENOME)

        assert not alias_dir.exists(), "the re-run left the alias tree behind"
        assert not build_tree.exists(), "the re-run left the build tree behind"
