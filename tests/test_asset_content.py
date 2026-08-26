"""
Write-path correctness for the asset content manager (refgenie.managers.asset).

``add_from_path`` is the only write path into the catalog -- build, pull, and
``Refgenie.add`` all route through it. Two families of invariant live here:

* **Interruption / write ordering.** Refgenie writes to three places that must
  agree -- the SQLite catalog, the RefgetStore, and plain files. Only the catalog
  has real transactions, so the rule is: commit the catalog first, then clean up
  the filesystem, and make the cleanup safe to re-run. Each test kills one write
  path between its filesystem phase and its SQL commit (or the reverse) and
  asserts the result is recoverable: no committed row references missing content;
  no content is destroyed by a rolled-back transaction; re-running succeeds with
  no ``--force`` and no hand-deleting. A ``SIGKILL`` runs no ``except`` block and
  no ``__exit__``, so where a path has a compensating handler the test disables it.

* **Rejected writes leave no trace.** ``data/`` has no sweeper, so content parked
  at a digest-addressed path with no catalog row is invisible, un-reapable, and
  satisfies the "destination already exists" branch of the next placement -- a
  later build silently adopts content nobody validated. A rejected ``add_from_path``
  must therefore leave the caller's source in place and place nothing.

The removal-focused ordering tests (asset delete before/after COMMIT,
genome-remove ordering) live in test_asset_removal.py.
"""

import os
import shutil
from pathlib import Path
from unittest.mock import MagicMock

import pytest
from sqlmodel import Session, select

from refgenie import Refgenie
from refgenie.db.tables import Asset, AssetName, StagedAsset
from refgenie.exceptions import MissingAssetError
from refgenie.managers.asset.colocation import (
    create_colocation_symlinks,
    get_colocation_filenames,
    get_colocation_metadata,
)
from refgenie.managers.asset.manager import AssetManager
from refgenie.utils.build import build_level1_to_digest
from refgenie.utils.staging import staged_archive_path, staged_archive_relpath
from tests.helpers import (
    ASSET,
    GENOME,
    GROUP,
    asset_name_rows,
    assets_rows,
    boom,
    build_rcrsd,
    fasta_asset,
    make_built_refgenie,
    make_engine,
    stage_copy,
)

# Component tier: real genome folders, real asset files, real tarballs. Deselected
# from the bare `pytest` inner loop; run with `pytest -m component`.
pytestmark = pytest.mark.component


def _group_dir(r: Refgenie, group: str) -> Path:
    """The data/<genome_digest>/<group> directory (for any group name)."""
    return r.genome_folder / "data" / r.alias.resolve(GENOME) / group


class TestBuildOrdering:
    """The completion flag must follow the commit, and re-runs must reclaim content."""

    def test_build_interrupted_before_commit_is_rebuildable(self, refgenie_fs, monkeypatch):
        """
        The completion flag used to be the last recipe command, so it appeared
        before ``add_from_path`` committed. A build killed in that gap left a flag
        asserting an asset that did not exist -- and both snakemake (declared output
        present) and pypiper (target exists) then skipped the rebuild forever.
        """
        r = refgenie_fs
        flag = r.asset._asset_builder._get_build_flag(GENOME, GROUP, ASSET)

        monkeypatch.setattr(r.asset, "add_from_path", boom)
        with pytest.raises(RuntimeError):
            r.build_asset(
                recipe_name="fasta",
                genome_name=GENOME,
                asset_group_name=GROUP,
                asset_name=ASSET,
            )

        assert not flag.exists(), "the flag was written for an asset that is not in the catalog"
        assert not r.asset.exists(
            genome_digest=r.alias.resolve(GENOME), asset_group_name=GROUP, asset_name=ASSET
        )

        monkeypatch.undo()
        asset = r.build_asset(
            recipe_name="fasta",
            genome_name=GENOME,
            asset_group_name=GROUP,
            asset_name=ASSET,
        )
        assert asset is not None, "the interrupted build could not be re-run"
        assert flag.exists(), "the flag was not written after the commit"

    def test_stale_build_flag_does_not_block_rebuild(self, refgenie_fs):
        """
        A flag with no catalog row behind it is debris, whoever left it there.
        Left in place it makes pypiper skip every command and hand an empty output
        folder to ``add_from_path``.
        """
        r = refgenie_fs
        flag = r.asset._asset_builder._get_build_flag(GENOME, GROUP, ASSET)
        flag.parent.mkdir(parents=True, exist_ok=True)
        flag.touch()

        asset = r.build_asset(
            recipe_name="fasta",
            genome_name=GENOME,
            asset_group_name=GROUP,
            asset_name=ASSET,
        )

        assert asset is not None, "a stale flag blocked the build"
        assert (r.genome_folder / asset.path / f"{r.alias.resolve(GENOME)}.fa").exists()

    def test_place_content_discards_source_destination_already_holds(self, tmp_path):
        """
        How a re-run of an interrupted build reclaims itself: the content moved to
        its digest path last time, the commit did not happen, so the rebuild
        produces byte-identical content, finds the destination occupied, and drops
        its own copy. Nothing is duplicated and nothing is lost.
        """
        source = tmp_path / "build_output"
        source.mkdir()
        (source / "f.txt").write_text("identical")
        destination = tmp_path / "data" / "digest"
        destination.mkdir(parents=True)
        (destination / "f.txt").write_text("identical")

        moved = AssetManager._place_content(source, destination)

        assert moved is False, "reported a move it did not make; the undo would delete real content"
        assert not source.exists(), "the redundant source was left behind"
        assert (destination / "f.txt").read_text() == "identical"


class TestStage:
    """Never publish a staged row for a directory that is not actually complete."""

    def test_stage_repairs_a_partially_linked_directory(self, refgenie_built):
        """
        File-mode staging creates the link directory and then links into it, one
        child at a time, and commits the row last. Killed mid-loop, a re-run used to
        see the directory, log "Skipping stage", and commit the row anyway --
        publishing a half-linked directory as complete.
        """
        r = refgenie_built
        asset = fasta_asset(r)
        content = r.genome_folder / asset.path
        expected = {p.name for p in content.iterdir()}

        link_dir = r.genome_stage_folder / r.alias.resolve(GENOME) / GROUP / ASSET
        link_dir.mkdir(parents=True)
        one = sorted(expected)[0]
        (link_dir / one).symlink_to(content / one)

        r.stage.create(asset, r.genome_folder, r.genome_stage_folder)

        linked = {p.name for p in link_dir.iterdir()}
        assert linked == expected, f"stage published an incomplete directory: {linked}"
        assert all((link_dir / name).resolve().exists() for name in linked)

    def test_stage_recreates_a_tarball_the_catalog_still_claims(self, refgenie_built):
        """
        A staged row is a claim that the tarball is there and is that digest, not
        evidence of it. Re-staging must check, or a missing tarball is published as
        a download that fails for every client.
        """
        r = refgenie_built
        asset = fasta_asset(r)
        staged = r.stage.create(asset, r.genome_folder, r.genome_stage_folder)
        archive = next(s for s in staged if s.mode == "archive")

        tarball = r.genome_stage_folder / r.alias.resolve(GENOME) / GROUP / f"{asset.digest}.tgz"
        assert tarball.is_file()
        original_digest = archive.tarball_digest
        tarball.unlink()

        restaged = r.stage.create(asset, r.genome_folder, r.genome_stage_folder)

        assert tarball.is_file(), "re-staging did not restore the missing tarball"
        repaired = next(s for s in restaged if s.mode == "archive")
        assert repaired.tarball_digest is not None
        assert repaired.tarball_size == tarball.stat().st_size
        # The archive is rebuilt from the same content, so the digest is stable.
        assert repaired.tarball_digest == original_digest


class TestRename:
    """A rename must purge both trees the old name owned, not just the alias tree."""

    def test_rename_purges_the_build_tree_of_the_old_name(self, refgenie_built):
        """
        ``rename`` cleared the stale alias tree but not ``builds/.../<old_name>/``.
        The completion flag left there asserts an asset that no longer answers to
        that name: a later build of the old name is skipped as already done, and a
        rename back to it fails with "already exists".
        """
        r = refgenie_built
        old_alias_dir = r.genome_folder / "alias" / GENOME / GROUP / ASSET
        old_build_dir = r.genome_folder / "builds" / GENOME / GROUP / ASSET
        assert old_alias_dir.exists() and old_build_dir.exists()

        r.asset.rename(
            asset_group_name=GROUP,
            asset_name=ASSET,
            new_asset_name="renamed",
            genome_name=GENOME,
        )

        assert not old_alias_dir.exists(), "the stale alias tree survived the rename"
        assert not old_build_dir.exists(), "the stale build tree survived the rename"
        assert (r.genome_folder / "alias" / GENOME / GROUP / "renamed").exists()


class TestPurge:
    """An interrupted purge must leave inert data, never a catalog that lies."""

    def test_purge_drops_the_schema_before_removing_the_files(self, refgenie_built, monkeypatch):
        """
        An interrupted purge should leave data with no catalog -- inert bytes -- not
        a catalog describing data that is already gone, which every read path would
        believe.
        """
        import refgenie.core.lifecycle as lifecycle_module

        r = refgenie_built
        # Read before the purge: afterwards there is no configuration table to read
        # the genome folder from -- which is the point.
        genome_folder = r.genome_folder
        monkeypatch.setattr(lifecycle_module.shutil, "rmtree", boom)

        with pytest.raises(RuntimeError):
            r.purge(force=True)

        assert not r.check_table_exists(), "the files were removed before the schema was dropped"
        assert (genome_folder / "data").exists(), "precondition: the files are still there"


class TestRejectedAddFromPath:
    """A rejected add_from_path must leave the caller's source in place and place nothing."""

    def test_cross_group_collision_leaves_the_source_in_place(self, refgenie_built):
        """
        Content whose digest is already claimed by another group is rejected.

        The rejection is decidable from the catalog alone, so it must happen before
        the content is moved: the destination is in a group the content does not
        belong to, and there is no row that would ever reference it.
        """
        r = refgenie_built
        staging, _ = stage_copy(r, fasta_asset(r), "staging_same")

        with pytest.raises(ValueError, match="Cross-group"):
            r.asset.add_from_path(
                asset_class_name="fasta",
                path=Path("staging_same"),
                asset_group_name="other_group",
                genome_name=GENOME,
                asset_name=ASSET,
            )

        assert staging.is_dir(), "rejected add consumed the caller's source directory"
        assert not _group_dir(r, "other_group").exists(), (
            f"rejected add left orphaned content: {list(_group_dir(r, 'other_group').rglob('*'))}"
        )

    def test_missing_seek_key_file_leaves_no_orphaned_content(self, refgenie_built):
        """
        Content that does not satisfy the asset class's declared seek keys is
        rejected while binding them -- after the placement, and after the row is
        built but before it is committed. The content must not survive at its
        digest path with no row pointing at it.
        """
        r = refgenie_built
        staging = r.genome_folder / "staging_incomplete"
        staging.mkdir()
        (staging / "not_a_fasta.txt").write_text("this satisfies no declared seek key\n")

        with pytest.raises(FileNotFoundError):
            r.asset.add_from_path(
                asset_class_name="fasta",
                path=Path("staging_incomplete"),
                asset_group_name="incomplete_group",
                genome_name=GENOME,
                asset_name="broken",
            )

        assert staging.is_dir(), "rejected add consumed the caller's source directory"
        assert not _group_dir(r, "incomplete_group").exists(), (
            f"rejected add left orphaned content: {list(_group_dir(r, 'incomplete_group').rglob('*'))}"
        )


# --- Pull re-runs: content with no row is debris, and debris is reclaimable ---
# These require fastapi and reuse the server/client harness from test_pull.py.


@pytest.fixture
def pull_setup(tmp_path, fixtures_path):
    """A real server app world (built + staged fasta asset) and a fresh client.

    Not the ``server_client_world`` fixture: these tests re-enter the serving
    context once per pull attempt.
    """
    from tests.helpers import make_server_client_world, requires_server

    requires_server()
    return make_server_client_world(make_engine(), tmp_path, fixtures_path)


def _pull(client_rg, server_rg, **kwargs):
    from tests.helpers import MOCK_SERVER_URL as url, serve_refgenie

    with serve_refgenie(client_rg, server_rg, url):
        return client_rg.pull(
            alias_name=GENOME,
            asset_group_name=GROUP,
            force_large=True,
            force_server_urls=[url],
            **kwargs,
        )


class TestPull:
    """A pull killed mid-write must be reclaimable on re-run, never wedged."""

    def test_pull_reclaims_an_orphaned_asset_directory(self, pull_setup, monkeypatch):
        """
        Killed after ``untar`` and before the commit, a pull left content at its
        digest path with no row. Every re-run then hit ``asset_dir.exists()`` and
        raised ``PullSkipped``, so only a human deleting the directory could clear
        it. ``PullTransaction`` is neutered here because a ``SIGKILL`` runs no
        ``__exit__``.
        """
        from refgenie.managers.asset.puller import PullTransaction

        server_rg, client_rg = pull_setup
        monkeypatch.setattr(PullTransaction, "_rollback", lambda self: None)
        monkeypatch.setattr(client_rg.asset, "add_from_path", boom)

        with pytest.raises(RuntimeError):
            _pull(client_rg, server_rg, force=False)

        genome_digest = client_rg.alias.resolve(GENOME)
        group_dir = client_rg.genome_folder / "data" / genome_digest / GROUP
        orphans = [p for p in group_dir.iterdir() if p.is_dir()]
        assert orphans, "precondition: the interrupted pull left content on disk"

        monkeypatch.undo()
        result = _pull(client_rg, server_rg, force=False)

        assert result is not None, "the orphan blocked every subsequent pull"
        assert client_rg.asset.exists(
            genome_digest=genome_digest, asset_group_name=GROUP, asset_name="default"
        )
        content = client_rg.genome_folder / result.path
        assert (content / f"{genome_digest}.fa").exists(), "the reclaimed pull has no content"

    def test_pull_re_renders_alias_tree_but_still_refuses_complete(self, pull_setup):
        """
        Killed after the commit and before ``_create_symlinks_for_alias``, a pull
        left the asset unreachable by name, and every re-run raised
        ``AssetExistsError`` before reaching any repair. The alias tree is derived
        from the catalog, so a re-pull renders it rather than refusing -- but a
        genuinely complete asset is still refused (the repair path is not a silent
        no-op success).
        """
        from refgenie.exceptions import AssetExistsError

        server_rg, client_rg = pull_setup

        assert _pull(client_rg, server_rg, force=False) is not None
        alias_dir = client_rg.genome_folder / "alias" / GENOME / GROUP / "default"
        assert alias_dir.exists()
        shutil.rmtree(alias_dir)

        # A missing (catalog-derived) alias tree is re-rendered, not refused.
        result = _pull(client_rg, server_rg, force=False)
        assert result is not None, "a missing alias tree made the pull unrepeatable"
        assert alias_dir.exists(), "the alias tree was not re-rendered"

        # But a genuinely complete asset is still refused.
        with pytest.raises(AssetExistsError):
            _pull(client_rg, server_rg, force=False)


# --- Colocation symlinks and staging ----------------------------------------


class TestColocation:
    """refgenie.managers.asset.colocation: parent-asset symlinks and metadata.

    Unit tier: mocks and tmp_path only, so it carries an explicit ``unit``
    mark against this module's ``component`` default.
    """

    pytestmark = pytest.mark.unit

    def _create_mock_fasta_asset(self, genome_folder: Path, genome_digest: str) -> Path:
        """Create a mock fasta asset directory with a .fa file."""
        fasta_dir = genome_folder / "data" / genome_digest / "fasta" / "default"
        fasta_dir.mkdir(parents=True)
        fasta_file = fasta_dir / f"{genome_digest}.fa"
        fasta_file.write_text(">chr1\nACGTACGT\n")
        return fasta_file

    def _make_mock_asset(self, seek_keys_dict: dict[str, Path]) -> MagicMock:
        """Create a mock Asset exposing a seek_keys_dict property."""
        asset = MagicMock()
        asset.seek_keys_dict = seek_keys_dict
        return asset

    def test_create_symlink_basic(self, tmp_path):
        """Colocation creates a *relative* symlink resolving to the parent file."""
        genome_folder = tmp_path / "genome_folder"
        genome_folder.mkdir()
        genome_digest = "abc123deadbeef"
        fasta_file = self._create_mock_fasta_asset(genome_folder, genome_digest)
        output_folder = genome_folder / "data" / genome_digest / "bwa_index" / "default"
        parent_asset = self._make_mock_asset(
            {"fasta": Path("data") / genome_digest / "fasta" / "default" / f"{genome_digest}.fa"}
        )
        input_assets = {
            "fasta": {
                "asset_class": "fasta",
                "default": "fasta",
                "colocate": [{"source_key": "fasta"}],
            }
        }

        created = create_colocation_symlinks(
            output_folder=output_folder,
            genome_folder=genome_folder,
            input_assets=input_assets,
            resolved_assets={"fasta": parent_asset},
        )

        assert created == [f"{genome_digest}.fa"]
        link_path = output_folder / f"{genome_digest}.fa"
        assert link_path.is_symlink()
        assert not os.path.isabs(os.readlink(link_path))  # relative link
        assert link_path.resolve() == fasta_file.resolve()

    def test_create_symlink_with_dest(self, tmp_path):
        """The `dest` field overrides the created symlink's name."""
        genome_folder = tmp_path / "genome_folder"
        genome_folder.mkdir()
        genome_digest = "abc123deadbeef"
        self._create_mock_fasta_asset(genome_folder, genome_digest)
        output_folder = genome_folder / "data" / genome_digest / "bwa_index" / "default"
        parent_asset = self._make_mock_asset(
            {"fasta": Path("data") / genome_digest / "fasta" / "default" / f"{genome_digest}.fa"}
        )
        input_assets = {
            "fasta": {
                "asset_class": "fasta",
                "default": "fasta",
                "colocate": [{"source_key": "fasta", "dest": "genome.fa"}],
            }
        }

        created = create_colocation_symlinks(
            output_folder=output_folder,
            genome_folder=genome_folder,
            input_assets=input_assets,
            resolved_assets={"fasta": parent_asset},
        )

        assert created == ["genome.fa"]
        link_path = output_folder / "genome.fa"
        assert link_path.is_symlink()
        assert link_path.resolve().is_file()

    @pytest.mark.parametrize(
        "input_assets,resolved",
        [
            ({"fasta": {"asset_class": "fasta", "default": "fasta"}}, {"fasta": MagicMock()}),
            (None, None),
        ],
    )
    def test_create_symlink_no_colocation(self, input_assets, resolved):
        """Recipes without colocate (or None inputs) create no symlinks."""
        created = create_colocation_symlinks(
            output_folder=Path("/tmp/test"),
            genome_folder=Path("/tmp"),
            input_assets=input_assets,
            resolved_assets=resolved,
        )
        assert created == []

    @pytest.mark.parametrize(
        "colocate,expected",
        [
            ([{"source_key": "fasta"}], ["abc.fa"]),
            ([{"source_key": "fasta", "dest": "genome.fa"}], ["genome.fa"]),
        ],
    )
    def test_get_colocation_filenames(self, colocate, expected):
        """Predicted name is `dest` if given, else the parent seek-key basename."""
        parent_asset = self._make_mock_asset({"fasta": Path("data/abc/fasta/default/abc.fa")})
        input_assets = {"fasta": {"asset_class": "fasta", "default": "fasta", "colocate": colocate}}
        assert get_colocation_filenames(input_assets, {"fasta": parent_asset}) == expected

    @pytest.mark.parametrize(
        "input_assets,expected",
        [
            (
                {
                    "fasta": {
                        "asset_class": "fasta",
                        "default": "fasta",
                        "colocate": [{"source_key": "fasta", "dest": "genome.fa"}],
                    }
                },
                [{"parent_asset_group": "fasta", "source_key": "fasta", "dest": "genome.fa"}],
            ),
            (
                {
                    "fasta": {
                        "asset_class": "fasta",
                        "default": "fasta",
                        "colocate": [{"source_key": "fasta"}],
                    }
                },
                [{"parent_asset_group": "fasta", "source_key": "fasta"}],
            ),
            ({"fasta": {"asset_class": "fasta"}}, None),
            (None, None),
        ],
    )
    def test_get_colocation_metadata(self, input_assets, expected):
        """Metadata reflects source_key + optional dest; None when no colocate."""
        assert get_colocation_metadata(input_assets) == expected


class TestStaging:
    """
    refgenie.utils.staging + StageManager. Component tier: exercises real
    genome folders, tarballs, and SQLite together (~0.2-0.5s apiece), so it is
    deselected from the bare `pytest` inner loop. Run with `pytest -m component`.
    """

    pytestmark = pytest.mark.component

    @pytest.mark.parametrize("stage_folder", [Path("/var/stage"), "/var/stage"])
    def test_staged_archive_path_helpers(self, stage_folder):
        """relpath and full path agree, accept str or Path, and never drift."""
        relpath = staged_archive_relpath("genome123", "fasta", "abc456")
        assert relpath == "genome123/fasta/abc456.tgz"
        full = staged_archive_path(stage_folder, "genome123", "fasta", "abc456")
        assert full == Path("/var/stage/genome123/fasta/abc456.tgz")
        assert full.relative_to(stage_folder) == Path(relpath)

    def test_stage_creation_and_removal(self, staged_refgenie):
        """Create → remove → the asset is gone and the list shrinks by its count."""
        r = staged_refgenie
        staged_assets = list(r.stage.list_all())
        assert len(staged_assets) >= 1
        staged, asset = staged_assets[0]
        assert isinstance(staged, StagedAsset)
        assert isinstance(asset, Asset)

        digest = staged.asset_digest
        same_digest_count = sum(
            1 for s, _ in staged_assets if s is not None and s.asset_digest == digest
        )

        r.stage.remove(asset_digest=digest)

        assert r.stage.exists_by_asset_digest(digest) is False
        remaining = list(r.stage.list_all())
        assert len(remaining) == len(staged_assets) - same_digest_count

    def test_manager_and_db_agree(self, staged_refgenie):
        """
        StageManager.list_all() count matches a direct DB query, and each staged
        row links to its asset by digest.
        """
        r = staged_refgenie
        with Session(r.database_engine) as session:
            db_staged = session.exec(select(StagedAsset)).all()

        manager_staged = list(r.stage.list_all())
        manager_count = len([s for s, _ in manager_staged if s is not None])

        assert len(db_staged) >= 1
        assert len(db_staged) == manager_count
        for staged, asset in manager_staged:
            if staged is not None:
                assert staged.asset_digest == asset.digest


# --- Digest-addressed asset names -------------------------------------------


OLD_NAME = "0.7.17"
NEW_NAME = "0.7.19"


@pytest.fixture
def refgenie_built_versioned(refgenie_fs):
    """Refgenie on a real filesystem with one FASTA asset built as ``0.7.17``.

    Deliberately NOT named ``refgenie_built``: conftest owns that name and
    builds the asset as ``test``, and shadowing it silently changed the asset
    name for every test in this module.
    """
    build_rcrsd(refgenie_fs, asset_name=OLD_NAME)
    return refgenie_fs


@pytest.fixture
def genome_digest(refgenie_built_versioned):
    return refgenie_built_versioned.alias.resolve("rCRSd")


class TestAssetNameContent:
    """Digest-addressed content with many names, and the alias-tree layout (component tier)."""

    pytestmark = pytest.mark.component

    def test_rebuild_identical_content_yields_one_asset_two_names(
        self, refgenie_built_versioned, genome_digest
    ):
        """A rebuild under a new name reuses the one content row and adds a name row;
        Asset.names exposes every name for use by the pull path."""
        r = refgenie_built_versioned
        build_rcrsd(r, asset_name=NEW_NAME)
        assert len(assets_rows(r)) == 1
        assert {n.name for n in asset_name_rows(r)} == {OLD_NAME, NEW_NAME}
        asset = r.asset.get(
            genome_digest=genome_digest, asset_group_name="fasta", asset_name=OLD_NAME
        )
        assert {n.name for n in asset.names} == {OLD_NAME, NEW_NAME}

    def test_both_names_resolve_across_the_api(self, refgenie_built_versioned, genome_digest):
        """Both names resolve via get/exists/seek/get_seek_key/set_parents to one digest."""
        r = refgenie_built_versioned
        build_rcrsd(r, asset_name=NEW_NAME)

        for name in (OLD_NAME, NEW_NAME):
            assert r.asset.exists(
                genome_digest=genome_digest, asset_group_name="fasta", asset_name=name
            )
            assert r.asset.get(
                genome_digest=genome_digest, asset_group_name="fasta", asset_name=name
            ) is not None
            assert r.asset.get_seek_key(
                genome_digest=genome_digest,
                asset_group_name="fasta",
                asset_name=name,
                seek_key_name="fasta",
            ) is not None
            assert r.asset.seek("rCRSd", "fasta", name, seek_key_name="fasta")

        a_old = r.asset.get(
            genome_digest=genome_digest, asset_group_name="fasta", asset_name=OLD_NAME
        )
        a_new = r.asset.get(
            genome_digest=genome_digest, asset_group_name="fasta", asset_name=NEW_NAME
        )
        assert a_old.digest == a_new.digest

        # set_parents resolves by either name without error.
        r.asset.set_parents(
            genome_digest=genome_digest,
            asset_group_name="fasta",
            asset_name=NEW_NAME,
            parent_asset_digests=[],
        )

    def test_content_directory_is_named_by_content_digest(self, refgenie_session):
        """The data directory is named by the content digest, not the asset name."""
        r = refgenie_session
        asset = r.asset.get(
            genome_digest=r.alias.resolve("rCRSd"), asset_group_name="fasta", asset_name="test"
        )
        assert Path(asset.path).name == asset.digest
        assert (r.genome_folder / asset.path).is_dir()

    def test_two_aliases_produce_full_cross_product_with_rewritten_filenames(
        self, refgenie_built_versioned, genome_digest
    ):
        """Adding a second alias renders every name x every alias, filenames rewritten."""
        r = refgenie_built_versioned
        build_rcrsd(r, asset_name=NEW_NAME)
        r.set_genome_alias(alias_name="rCRSd_v2", genome_digest=genome_digest)

        for alias in ("rCRSd", "rCRSd_v2"):
            for name in (OLD_NAME, NEW_NAME):
                d = r.alias_folder / alias / "fasta" / name
                assert d.is_dir(), f"missing alias dir {d}"
                # Filenames are rewritten to the owning alias (regression guard).
                fa = d / f"{alias}.fa"
                assert fa.exists(), f"missing rewritten file {fa}"
                # The digest-named file must NOT appear under the alias tree.
                assert not (d / f"{genome_digest}.fa").exists()

        # Every name is a view onto the same underlying content file.
        resolved = {
            (r.alias_folder / "rCRSd" / "fasta" / name / "rCRSd.fa").resolve()
            for name in (OLD_NAME, NEW_NAME)
        }
        assert len(resolved) == 1
        assert resolved.pop().is_file()

    def test_seek_returns_alias_path_and_abs_returns_content(self, refgenie_session):
        """seek returns an existing alias path; abs_path=True returns the content path."""
        r = refgenie_session
        alias_path = r.asset.seek("rCRSd", "fasta", "test", seek_key_name="fasta", force_exists=True)
        content_path = r.asset.seek(
            "rCRSd", "fasta", "test", seek_key_name="fasta", force_exists=True, abs_path=True
        )
        assert "/alias/" in alias_path
        assert "/data/" in content_path
        assert alias_path != content_path
        assert Path(alias_path).read_bytes() == Path(content_path).read_bytes()

    def test_supplied_digest_becomes_the_asset_identity(self, refgenie_built_versioned, genome_digest):
        """A caller-supplied digest keys the asset verbatim rather than being recomputed.

        This is the pull path: which files the digest covers is decided by the
        building recipe's `inherent` set, which does not travel with the archive,
        so a client that recomputed would disagree with the server.
        """
        r = refgenie_built_versioned
        asset = r.asset.get(
            genome_digest=genome_digest, asset_group_name="fasta", asset_name=OLD_NAME
        )
        _, staging = stage_copy(r, asset, "staging_supplied")

        supplied = "f" * 64
        added = r.asset.add_from_path(
            asset_class_name="fasta",
            path=staging,
            asset_group_name="fasta",
            genome_name="rCRSd",
            asset_name="pulled",
            digest=supplied,
        )

        assert added.digest == supplied, "the supplied digest must be used verbatim"
        assert r.asset.get_by_digest(supplied) is not None


class TestBuildProvenanceOnNames:
    """Build provenance lives on the name row, not the content row (component tier).

    Before this, a build's provenance was written as seek keys on the ``asset``
    row. A second build producing byte-identical content takes the reconcile
    path, which adds only an ``AssetName`` row and never touches the content's
    seek keys -- so the second build's provenance was silently dropped. These
    tests pin the invariant that replaced it: one row per build, on the name.
    """

    pytestmark = pytest.mark.component

    def test_a_build_records_its_provenance_on_the_name_row(self, refgenie_built_versioned):
        r = refgenie_built_versioned
        (row,) = asset_name_rows(r)

        assert row.build_digest is not None
        assert row.build_timestamp is not None
        assert row.refgenie_version
        assert row.recipe_id is not None
        # Re-derivable from the row alone, with no access to the inputs.
        assert build_level1_to_digest(row.build_level1) == row.build_digest

    def test_second_build_of_identical_content_keeps_the_first_provenance(
        self, refgenie_built_versioned
    ):
        """The regression test for the bug the provenance columns fix.

        The rebuild is the same recipe, genome, version and params, so it is
        the same *build* -- it gets its own name and its own full provenance,
        without disturbing the first name's row.
        """
        r = refgenie_built_versioned
        (first,) = asset_name_rows(r)
        original = (first.build_digest, first.build_timestamp, first.refgenie_version)

        build_rcrsd(r, asset_name=NEW_NAME)

        rows = {row.name: row for row in asset_name_rows(r)}
        assert set(rows) == {OLD_NAME, NEW_NAME}
        assert len(assets_rows(r)) == 1, "identical content must stay one asset row"
        assert (
            rows[OLD_NAME].build_digest,
            rows[OLD_NAME].build_timestamp,
            rows[OLD_NAME].refgenie_version,
        ) == original, "the first build's provenance must survive the second build"
        # Same build, two names: both record it, since unique_group_build_digest
        # is a non-unique index precisely so this does not collide.
        assert rows[NEW_NAME].build_digest == rows[OLD_NAME].build_digest
        assert rows[NEW_NAME].build_timestamp is not None
        assert rows[NEW_NAME].build_timestamp != rows[OLD_NAME].build_timestamp

    def test_two_genomes_build_fasta_without_colliding(self, tmp_path, fixtures_path):
        """The regression test for a digest that omits the genome.

        ``fasta_asset_recipe.yaml`` declares no input files, params or assets, so
        every genome's fasta build would compute one digest and the second
        insert would die on the unique index.
        """
        r = make_built_refgenie(tmp_path / "two_genomes", fixtures_path, build=True)
        r.genome.initialize_genome(
            fasta_file_path=fixtures_path / "demo.fa",
            alias_names=["demo"],
            description="demo genome",
        )
        r.build_asset(
            recipe_name="fasta",
            genome_name="demo",
            asset_group_name=GROUP,
            asset_name=ASSET,
        )

        digests = [row.build_digest for row in asset_name_rows(r) if row.build_digest]
        assert len(digests) == 2, "both genomes must record a build"
        assert len(set(digests)) == 2, "two genomes' fasta builds must not share a digest"


class TestAssetNameGuards:
    """Name/content collision and removal guards (component tier)."""

    pytestmark = pytest.mark.component

    def test_claiming_existing_name_for_different_content_raises(
        self, refgenie_built_versioned, genome_digest
    ):
        """A name already mapped to content cannot be reused for different content."""
        r = refgenie_built_versioned
        asset = r.asset.get(
            genome_digest=genome_digest, asset_group_name="fasta", asset_name=OLD_NAME
        )
        _, staging = stage_copy(r, asset, "staging_diff", mutate=b"\nN\n")

        with pytest.raises(ValueError):
            r.asset.add_from_path(
                asset_class_name="fasta",
                path=staging,
                asset_group_name="fasta",
                genome_name="rCRSd",
                asset_name=OLD_NAME,
            )

    def test_remove_via_noncanonical_name_raises(self, refgenie_built_versioned, genome_digest):
        """Removal via a secondary name raises rather than silently dropping siblings."""
        r = refgenie_built_versioned
        build_rcrsd(r, asset_name=NEW_NAME)
        canonical = r.asset.get(
            genome_digest=genome_digest, asset_group_name="fasta", asset_name=OLD_NAME
        ).name
        secondary = NEW_NAME if canonical == OLD_NAME else OLD_NAME

        with pytest.raises(ValueError):
            r.asset.remove(
                genome_digest=genome_digest, asset_group_name="fasta", asset_name=secondary
            )


class TestAssetNameDefaults:
    """Default-flag invariants across build and add_from_path (component tier)."""

    pytestmark = pytest.mark.component

    def test_set_default_records_exact_name_and_keeps_one_flag(
        self, refgenie_built_versioned, genome_digest
    ):
        """set_default records the name the caller passed (not the canonical name) and
        flipping it clears the previous flag -- at most one default per group."""
        r = refgenie_built_versioned
        build_rcrsd(r, asset_name=NEW_NAME)
        for name in (NEW_NAME, OLD_NAME):
            r.asset.set_default(
                genome_digest=genome_digest, asset_group_name="fasta", asset_name=name
            )
            assert r.asset.get_default("fasta", genome_digest=genome_digest) == name
            with Session(r.database_engine) as session:
                defaults = session.exec(select(AssetName).where(AssetName.is_default)).all()
            assert [d.name for d in defaults] == [name]

    def test_set_default_from_different_group_raises(self, refgenie_built_versioned, genome_digest):
        """A name that belongs to no group cannot be this group's default."""
        r = refgenie_built_versioned
        with pytest.raises(MissingAssetError):
            r.asset.set_default(
                genome_digest=genome_digest,
                asset_group_name="fasta",
                asset_name="a-name-from-nowhere",
            )

    def test_build_promotes_new_version_to_default(self, refgenie_built_versioned, genome_digest):
        """A deliberate build of a new version becomes the group default."""
        r = refgenie_built_versioned
        assert r.asset.get_default("fasta", genome_digest=genome_digest) == OLD_NAME
        build_rcrsd(r, asset_name=NEW_NAME)
        assert r.asset.get_default("fasta", genome_digest=genome_digest) == NEW_NAME

    @pytest.mark.parametrize(
        "set_default, group, expected",
        [
            # The pull path (None) must not hijack an existing default.
            pytest.param(None, "fasta", OLD_NAME, id="none-preserves-existing"),
            pytest.param(True, "fasta", "pulled", id="true-promotes"),
            # False never promotes, even when the group is brand new (where
            # None's legacy behavior would have promoted).
            pytest.param(False, "fasta_new", None, id="false-never-promotes"),
        ],
    )
    def test_add_from_path_set_default(
        self, refgenie_built_versioned, genome_digest, set_default, group, expected
    ):
        r = refgenie_built_versioned
        asset = r.asset.get(
            genome_digest=genome_digest, asset_group_name="fasta", asset_name=OLD_NAME
        )
        _, staging = stage_copy(r, asset, "staging_default", mutate=b"\nN\n")

        r.asset.add_from_path(
            asset_class_name="fasta",
            path=staging,
            asset_group_name=group,
            genome_name="rCRSd",
            asset_name="pulled",
            set_default=set_default,
        )

        assert r.asset.get_default(group, genome_digest=genome_digest) == expected
