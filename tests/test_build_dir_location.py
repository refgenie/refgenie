"""
Tests for the location and lifetime of the build bookkeeping directory.

Build bookkeeping (log, commands, profile, stats, and the completion flag) lives
in a top-level ``builds/`` tree keyed by build invocation, NOT inside the asset
directory it produced. Two properties follow, and both are load-bearing:

1. **Separation.** Nothing in the asset directory is build bookkeeping, so
   digests, tarballs, rsyncs, ``directory_contents``, and file-mode staging need
   no exclusion rules. This is enforced by construction rather than by five call
   sites each remembering to exclude a magic directory name.

2. **Lifetime.** The completion flag is not a historical record — it asserts
   that the asset currently exists. Removal must destroy it. If it survives,
   snakemake sees its declared output present, skips the rebuild, and the
   removed asset silently never comes back.
"""

import pytest
from sqlmodel import Session, select

from refgenie.const import ALIAS_DIR, BUILDS_DIR
from refgenie.db.tables import AssetClass
from refgenie.exceptions import MissingBuildInputError
from refgenie.utils.build import get_build_dir
from tests.helpers import (
    ASSET,
    GENOME,
    GROUP,
    build_dir,
    build_flag,
    content_dir,
    make_command_values,
)


# Component tier: these exercise real genome folders, real tarballs and real
# SQLite together, at ~0.2-0.5s apiece. Deselected from the bare `pytest`
# inner loop; run with `pytest -m component`.
pytestmark = pytest.mark.component


# --- Layout -----------------------------------------------------------------


def test_build_dir_is_outside_the_asset_directory(refgenie_built):
    """Build bookkeeping is written to builds/, never into the asset dir."""
    r = refgenie_built
    bdir = build_dir(r)
    assert bdir.is_dir(), f"build directory missing: {bdir}"
    assert bdir.is_relative_to(r.genome_folder / BUILDS_DIR)

    asset_dir = content_dir(r)
    assert not bdir.is_relative_to(asset_dir)
    # No build artifacts of any kind reached the asset directory.
    leaked = [
        p.name
        for p in asset_dir.rglob("*")
        if p.name.endswith((".flag", "_commands.sh", "_profile.tsv", "log.md"))
    ]
    assert leaked == [], f"build artifacts leaked into the asset directory: {leaked}"


def test_advertised_flag_is_the_real_flag(refgenie_built):
    """
    get_asset_build_target_template() points at the file build() actually wrote.

    Previously the template advertised an alias path while build() wrote a
    digest path, and a symlink bridged the two. Making them one path is the
    point of the relocation — a symlink here means the roots can drift again.
    """
    r = refgenie_built
    flag = build_flag(r)
    assert flag.is_file(), f"build flag missing: {flag}"
    assert not flag.is_symlink(), f"build flag should be a real file: {flag}"
    assert flag.parent == build_dir(r)


def test_file_mode_staging_excludes_build_artifacts(refgenie_built, fixtures_path):
    """
    File-mode staging symlinks every child of the asset dir into the served
    stage dir. With bookkeeping relocated there is nothing there to leak, so
    build logs and commands can no longer be published to a served bucket.
    """
    r = refgenie_built
    with Session(r.database_engine) as session:
        ac = session.exec(select(AssetClass).where(AssetClass.name == GROUP)).first()
        ac.serving_modes = ["file"]
        session.add(ac)
        session.commit()

    asset = list(r.asset.list_assets())[0]
    r.stage.create(
        asset=asset,
        genome_folder=r.genome_folder,
        genome_stage_folder=r.genome_stage_folder,
        build_dir=build_dir(r),
    )

    link_dir = (
        r.genome_stage_folder
        / asset.asset_group.genome_digest
        / asset.asset_group.name
        / asset.name
    )
    assert link_dir.is_dir()
    staged_names = [p.name for p in link_dir.iterdir()]
    assert not any(
        n.endswith((".flag", "_commands.sh", "_profile.tsv", "log.md")) for n in staged_names
    ), f"build artifacts staged for serving: {staged_names}"
    assert BUILDS_DIR not in staged_names
    assert "_refgenie_build" not in staged_names


def test_stage_records_build_commands_from_the_build_dir(refgenie_built):
    """The build commands are still captured, now read from builds/."""
    r = refgenie_built
    asset = list(r.asset.list_assets())[0]
    staged = r.stage.create(
        asset=asset,
        genome_folder=r.genome_folder,
        genome_stage_folder=r.genome_stage_folder,
        build_dir=build_dir(r),
    )
    assert staged
    assert staged[0].build_commands, "build commands were not recorded"


# --- Lifetime ---------------------------------------------------------------
# The flag-destruction and group-directory-cleanup contracts of asset removal
# now live in test_asset_removal.py, alongside the other removal tests.


def test_alias_delete_handler_sweeps_its_build_subtree(refgenie_built):
    """
    The builds/ tree is alias-keyed, so deleting an alias must sweep its whole
    subtree — the alias directory and the build directory alike.

    NOTE: this exercises the SQL-backed alias path (``AliasManager``, whose
    deletes fire ``_alias_before_delete_handler``). The default local mode uses
    ``StoreAliasManager``, which performs no filesystem cleanup at all — under
    that manager ``alias/<name>/`` is already orphaned on removal today, and
    ``builds/<name>/`` now orphans with it. That gap is pre-existing and lives
    in the store alias manager, not here.
    """
    from refgenie.db.events import _alias_before_delete_handler
    from refgenie.db.tables import Alias

    r = refgenie_built
    alias_build_dir = r.genome_folder / BUILDS_DIR / GENOME
    alias_dir = r.alias_folder / GENOME
    assert alias_build_dir.is_dir()
    assert alias_dir.is_dir()

    with Session(r.database_engine) as session:
        _alias_before_delete_handler(None, session.connection(), Alias(name=GENOME))

    assert not alias_build_dir.exists(), f"orphaned build subtree: {alias_build_dir}"
    assert not alias_dir.exists(), f"orphaned alias subtree: {alias_dir}"


def test_rebuild_after_removal_recreates_the_flag(refgenie_built, fixtures_path):
    """
    End-to-end: remove, rebuild, and the flag is back at the advertised path.

    This is the behavior snakemake depends on — a removed asset must look
    un-built, and rebuilding must satisfy the declared output.
    """
    r = refgenie_built
    flag = build_flag(r)

    r.asset.remove(
        genome_name=GENOME,
        asset_group_name=GROUP,
        asset_name=ASSET,
        keep_asset_group=True,
    )
    assert not flag.exists()

    r.build_asset(
        recipe_name="fasta",
        genome_name=GENOME,
        asset_group_name=GROUP,
        asset_name=ASSET,
    )
    assert flag.is_file(), f"rebuild did not recreate the flag: {flag}"


def test_rebuilding_an_existing_asset_leaves_the_flag_in_place(refgenie_built):
    """
    Re-running a build for an already-built asset skips the build and returns
    early. The flag must already be present, or snakemake raises
    MissingOutputException — the case the deleted symlink code existed to patch.
    """
    r = refgenie_built
    flag = build_flag(r)
    assert flag.is_file()

    r.build_asset(
        recipe_name="fasta",
        genome_name=GENOME,
        asset_group_name=GROUP,
        asset_name=ASSET,
    )

    assert flag.is_file(), "flag missing after a skipped rebuild"
    assert not flag.is_symlink()


def test_skipped_build_restores_a_missing_flag(refgenie_built):
    """
    builds/ has no sweeper, so it can be pruned independently of the assets it
    describes. If that happens, a rebuild takes the "already exists -> skip"
    early return and must still restore the declared output — otherwise every
    previously-built asset fails the next snakemake run with
    MissingOutputException.
    """
    r = refgenie_built
    flag = build_flag(r)
    flag.unlink()
    assert not flag.exists()

    r.build_asset(
        recipe_name="fasta",
        genome_name=GENOME,
        asset_group_name=GROUP,
        asset_name=ASSET,
    )

    assert flag.is_file(), "skipped build did not restore the missing flag"


# --- Alias-keying is not the alias the caller happens to name ---------------


def test_staging_finds_the_build_dir_under_a_different_alias(refgenie_built):
    """
    builds/ is keyed by the alias used at BUILD time. A genome can have several
    aliases, and a later `stage` may name a different one. Deriving the build
    path from the alias the user typed silently misses the bookkeeping and
    stages the asset with empty build_commands, losing provenance with no error.
    """
    r = refgenie_built
    genome_digest = r.alias.resolve(GENOME)
    r.alias.add("GRCh_other", genome_digest)

    # The asset was built as GENOME; look it up via the other alias.
    found = r.find_build_dir(
        genome_digest=genome_digest,
        asset_group_name=GROUP,
        asset_name=ASSET,
    )
    assert found == build_dir(r), "build dir not found via a non-build alias"


def test_find_build_dir_returns_none_for_an_asset_with_no_build_dir(refgenie_built):
    """A pulled asset has no build bookkeeping; staging must degrade quietly."""
    r = refgenie_built
    assert (
        r.find_build_dir(
            genome_digest=r.alias.resolve(GENOME),
            asset_group_name=GROUP,
            asset_name="never_built",
        )
        is None
    )


# --- Alias removal in local (store-backed) mode -----------------------------


def test_store_alias_removal_sweeps_the_build_tree(refgenie_built):
    """
    Local mode is the default and uses StoreAliasManager, which stores aliases
    in the refget store rather than the database — so the Alias before_delete
    event never fires and cannot do the cleanup. Without an explicit sweep,
    every alias removal in the default configuration orphans both alias/<name>/
    and builds/<name>/.
    """
    r = refgenie_built
    alias_dir = r.genome_folder / "alias" / GENOME
    genome_build_dir = get_build_dir(genome_folder=r.genome_folder, genome_name=GENOME)
    assert alias_dir.is_dir()
    assert genome_build_dir.is_dir()

    r.alias.remove(GENOME)

    assert not alias_dir.exists(), "alias directory orphaned by alias removal"
    assert not genome_build_dir.exists(), "build directory orphaned by alias removal"


def test_genome_removal_leaves_no_orphans_or_dangling_symlinks(refgenie_built):
    """
    Genome removal delegates to the alias manager per alias (GenomeManager.remove),
    so it inherits whatever that manager does — which is why the local-mode gap
    reached it too.

    This case is the worst of the two: genome removal DOES delete the data
    directory, so an un-swept alias tree is left full of symlinks pointing at
    files that no longer exist. Assert on dangling links directly, not just on
    directory existence.
    """
    r = refgenie_built
    genome_digest = r.alias.resolve(GENOME)
    alias_dir = r.genome_folder / ALIAS_DIR / GENOME
    genome_build_dir = get_build_dir(genome_folder=r.genome_folder, genome_name=GENOME)
    assert alias_dir.is_dir()
    assert genome_build_dir.is_dir()

    r.genome.remove(genome_digest)

    assert not alias_dir.exists(), "alias directory orphaned by genome removal"
    assert not genome_build_dir.exists(), "build directory orphaned by genome removal"
    dangling = [p for p in r.genome_folder.rglob("*") if p.is_symlink() and not p.exists()]
    assert dangling == [], f"genome removal left dangling symlinks: {dangling}"


def test_store_alias_manager_without_a_genome_folder_does_not_raise():
    """
    genome_folder_getter is optional on StoreAliasManager. One constructed
    without it (embedding, tests) must still remove the alias from the store
    rather than crashing on the cleanup it cannot perform.
    """
    from refgenie.managers.alias import StoreAliasManager

    removed: list[tuple[str, str]] = []

    class _FakeStore:
        def remove_collection_alias(self, namespace, name):
            removed.append((namespace, name))

        def get_collection_metadata_by_alias(self, namespace, name):
            return None

    mgr = StoreAliasManager(refget_store_getter=lambda: _FakeStore())
    mgr._cache["rCRSd"] = "somedigest"

    mgr.remove("rCRSd")

    assert removed == [("refgenie", "rCRSd")]
    assert "rCRSd" not in mgr._cache


# --- Pre-build input validation (unit tier) ----------------------------------


class FakeRecipe:
    """Minimal recipe stand-in for ``AssetBuilder._validate_build_inputs``."""

    def __init__(self, input_files=None, input_params=None):
        self.input_files = input_files
        self.input_params = input_params


class TestBuildInputValidation:
    """Pre-build input validation catches missing inputs (unit tier: fakes + tmp_path).

    Explicitly ``unit`` against this module's ``component`` default: it calls a
    static method with fakes and tmp_path and builds nothing.
    """

    pytestmark = pytest.mark.unit

    def test_missing_refget_store_detected(self, tmp_path):
        """Build validation catches a missing RefgetStore path."""
        from refgenie.managers.asset.builder import AssetBuilder

        command_values = make_command_values(
            genome_digest="fake_digest",
            asset_group_name="fasta",
            genome_folder=tmp_path / "nonexistent",
            refget_store_path=str(tmp_path / "nonexistent" / ".refget_store"),
        )

        with pytest.raises(MissingBuildInputError, match="RefgetStore not found"):
            AssetBuilder._validate_build_inputs(FakeRecipe(), command_values)

    def test_missing_required_file_detected(self, tmp_path):
        """Build validation catches missing required input files."""
        from refgenie.managers.asset.builder import AssetBuilder

        command_values = make_command_values(
            genome_digest="fake_digest",
            asset_group_name="bowtie2_index",
            genome_folder=tmp_path,
            refget_store_path=str(tmp_path / ".refget_store"),
        )
        (tmp_path / ".refget_store").mkdir(parents=True, exist_ok=True)

        recipe = FakeRecipe(input_files={"fasta": {"description": "Input FASTA file"}})
        with pytest.raises(MissingBuildInputError, match="Missing required file 'fasta'"):
            AssetBuilder._validate_build_inputs(recipe, command_values)

    def test_file_with_default_not_flagged(self, tmp_path):
        """Build validation does not flag input files that have defaults."""
        from refgenie.managers.asset.builder import AssetBuilder

        command_values = make_command_values(
            genome_digest="fake_digest",
            asset_group_name="test",
            genome_folder=tmp_path,
            refget_store_path=str(tmp_path / ".refget_store"),
        )
        (tmp_path / ".refget_store").mkdir(parents=True, exist_ok=True)

        recipe = FakeRecipe(
            input_files={"optional_file": {"description": "Optional", "default": "/some/path"}}
        )
        # Should not raise.
        AssetBuilder._validate_build_inputs(recipe, command_values)
