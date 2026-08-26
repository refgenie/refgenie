"""
Tests for the genome manager (refgenie.managers.genome / sources.genomes) and the
genome init/build flow.

Unit-tier classes cover catalog queries, getseq coordinate semantics, removal
cascade, build-input validation (pure static-method checks), and the genome-init
failure/exit contract (via lightweight doubles). Component-tier classes exercise
the real initialize_and_build flow, alias self-repair on re-init, and the
metadata-only store regression on a real on-disk gtars store.

Also here: Refgenie.init() folder-creation semantics, the alias manager and
its symlink tree, and getseq() against a RefgetStore backend. The ORM get
round-trip stays in test_db.py; it is not duplicated here.
"""

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest
from sqlmodel import select

from refgenie import Refgenie
from refgenie.cli.commands.genome import GenomeInitModel, handle_genome_init
from refgenie.core.sequences import _parse_locus
from refgenie.db.tables import Configuration
from refgenie.exceptions import MissingAliasError, MissingGenomeError
from refgenie.utils.symlinks import create_alias_symlinks
from tests.helpers import GENOME, add_asset_from_files, genome_digest, make_engine


class TestGenomeCatalog:
    """Catalog-level genome queries (unit tier)."""

    def test_list_all_and_exists(self, refgenie_session):
        """The built genome is listed and exists() agrees; unknown digests do not."""
        digest = refgenie_session.alias.resolve("rCRSd")
        assert digest in {g.digest for g in refgenie_session.genome.list_all()}
        assert refgenie_session.genome.exists(digest)
        assert not refgenie_session.genome.exists("nonexistent_genome")


class TestGenomeSequence:
    """getseq coordinate semantics and error paths (unit tier)."""

    def test_getseq_missing_chromosome_raises(self, refgenie_session):
        """A locus on a chromosome absent from the genome raises.

        Coordinate parsing and substring semantics are covered in the
        RefgetStore getseq section below; this missing-chromosome path is
        unique here.
        """
        with pytest.raises(ValueError):
            refgenie_session.getseq("rCRSd", "chr999:0-10")


class TestGenomeRemove:
    """Genome removal cascades to aliases; removing a missing genome raises (unit tier)."""

    def test_remove_cascades_to_aliases(self, refgenie_minimal):
        digest = "eeeeeeee1234567890abcdef12345678"
        aliases = ["remove_test_alias1", "remove_test_alias2"]
        refgenie_minimal.genome.add(
            digest=digest, description="Test genome for removal", alias_names=aliases
        )
        assert refgenie_minimal.genome.exists(digest)
        for alias in aliases:
            assert refgenie_minimal.alias.resolve(alias) == digest

        refgenie_minimal.genome.remove(digest)

        assert not refgenie_minimal.genome.exists(digest)
        for alias in aliases:
            with pytest.raises(MissingAliasError):
                refgenie_minimal.alias.resolve(alias)

    def test_remove_nonexistent_raises(self, refgenie_minimal):
        with pytest.raises(MissingGenomeError):
            refgenie_minimal.genome.remove("nonexistent_digest_12345678")


# --- genome-init failure/exit contract: doubles + tests (unit tier) ---------
#
# Regression: `build_asset` signals failure by RETURNING None rather than
# raising, so `handle_genome_init` wrapping the call in
# try/except caught nothing -- a failed build logged "built successfully" and
# exited 0. The orchestrator gates on that exit status, so it recorded the genome
# as initialized and the real error surfaced later as a misleading "not found".


def _init_cmd(**overrides):
    """Stand-in for the parsed `genome init --fasta ... --name ...` command."""
    fields = dict(
        name=["athaliana"], fasta="/tmp/does-not-matter.fa", server=None, store=None,
        namespace=None, digest=None, description="", species=None, fhr=None,
        force=False, build=True,
    )
    fields.update(overrides)
    return SimpleNamespace(**fields)


# The double must not drift from the real parsed model.
assert set(_init_cmd().__dict__) <= set(GenomeInitModel.model_fields)


def _init_refgenie(build_result):
    """A refgenie double whose build_asset outcome the caller controls."""
    rg = MagicMock()
    rg.genome.initialize_genome.return_value = ("AHyVecap0miRdrz6ZffizkZQUOJJBO2c", True)
    rg.recipe.exists.return_value = True
    if isinstance(build_result, Exception):
        rg.build_asset.side_effect = build_result
    else:
        rg.build_asset.return_value = build_result
    return rg


class TestGenomeInitFailure:
    """genome init must not report success when its automatic fasta build fails."""

    @pytest.mark.parametrize(
        "build_result",
        [None, RuntimeError("Collection not found")],
        ids=["build-returns-none", "build-raises"],
    )
    def test_failed_build_aborts_init(self, caplog, build_result):
        """A failed fasta build must exit non-zero and must not claim success."""
        refgenie = _init_refgenie(build_result)
        with caplog.at_level("INFO"):
            with pytest.raises(SystemExit) as excinfo:
                handle_genome_init(_init_cmd(), refgenie)
        assert excinfo.value.code == 1
        refgenie.build_asset.assert_called_once()
        assert "built successfully" not in caplog.text

    def test_succeeds_when_build_returns_asset(self):
        """A genuine success must NOT exit; that is what makes the failure cases meaningful."""
        refgenie = _init_refgenie(object())
        handle_genome_init(_init_cmd(), refgenie)
        refgenie.build_asset.assert_called_once()


@pytest.fixture
def init_t7(refgenie_with_fasta, fixtures_path):
    """Run initialize_and_build for the 't7' genome from the rCRSd fixture FASTA."""
    def _init(**overrides):
        kwargs = dict(
            fasta_file_path=fixtures_path / "rCRSd.fa",
            genome_names=["t7"],
            description="Test genome",
        )
        kwargs.update(overrides)
        return refgenie_with_fasta.initialize_and_build(**kwargs)
    return _init


@pytest.fixture
def reinit_t7(refgenie_with_fasta, fixtures_path):
    """Re-run initialize_genome for 't7' with use_existing=True (what --force passes)."""
    def _reinit():
        return refgenie_with_fasta.genome.initialize_genome(
            fasta_file_path=fixtures_path / "rCRSd.fa",
            description="Test genome",
            alias_names=["t7"],
            use_existing=True,
        )
    return _reinit


class TestInitializeAndBuild:
    """The initialize_and_build convenience method (component tier)."""

    pytestmark = pytest.mark.component

    def test_init_with_fasta_auto_builds(self, refgenie_with_fasta, init_t7):
        """initialize_and_build creates both the genome and the fasta asset."""
        digest, created = init_t7()
        assert created is True
        assert digest is not None

        genome_digest = refgenie_with_fasta.alias.resolve("t7")
        asset = refgenie_with_fasta.asset.get(
            genome_digest=genome_digest,
            asset_group_name="fasta",
            asset_name="default",
        )
        assert asset is not None

    def test_init_with_no_build(self, refgenie_with_fasta, init_t7):
        """build_fasta=False skips the fasta asset; it can be built separately after."""
        digest, created = init_t7(build_fasta=False)
        assert created is True
        assert digest is not None

        genome_digest = refgenie_with_fasta.alias.resolve("t7")
        from refgenie.exceptions import MissingAssetError

        with pytest.raises(MissingAssetError):
            refgenie_with_fasta.asset.get(
                genome_digest=genome_digest,
                asset_group_name="fasta",
                asset_name="default",
            )

        # Building fasta separately after a no-build init succeeds.
        refgenie_with_fasta.build_asset(
            recipe_name="fasta",
            genome_name="t7",
            asset_group_name="fasta",
        )
        asset = refgenie_with_fasta.asset.get(
            genome_digest=genome_digest,
            asset_group_name="fasta",
            asset_name="default",
        )
        assert asset is not None

    def test_init_idempotency_with_existing(self, init_t7):
        """Calling initialize_and_build twice with use_existing=True does not crash."""
        digest1, created1 = init_t7()
        assert created1 is True

        digest2, _created2 = init_t7(use_existing=True)
        assert digest1 == digest2


class TestAliasRepairOnReinit:
    """`genome init --force` must restore an alias that went missing (component tier).

    Genome rows live in the SQL catalog; aliases live in the RefgetStore, whose
    index is rewritten wholesale with no locking. A concurrent writer can drop an
    alias while leaving the genome row intact. Re-initialization short-circuits
    whenever the DIGEST row exists, so the retry must still repair the alias.
    """

    pytestmark = pytest.mark.component

    def test_reinit_restores_dropped_alias(self, refgenie_with_fasta, init_t7, reinit_t7):
        digest, created = init_t7(build_fasta=False)
        assert created is True
        assert refgenie_with_fasta.alias.resolve("t7") == digest

        # Simulate the race: alias gone, genome row still present.
        refgenie_with_fasta.alias.remove("t7")
        with pytest.raises(MissingAliasError):
            refgenie_with_fasta.alias.resolve("t7")
        assert refgenie_with_fasta.genome.exists(digest)

        # Re-init with use_existing (what `--force` passes) must repair it.
        digest2, created2 = reinit_t7()

        assert digest2 == digest
        assert created2 is False, "genome row already existed, so nothing was created"
        assert refgenie_with_fasta.alias.resolve("t7") == digest, (
            "re-init left the alias unresolvable, so the genome can never self-heal"
        )

    def test_reinit_leaves_existing_alias_alone(self, refgenie_with_fasta, init_t7, reinit_t7):
        """The repair must be a no-op when the alias is already correct."""
        digest, _ = init_t7(build_fasta=False)

        digest2, created2 = reinit_t7()

        assert digest2 == digest
        assert created2 is False
        assert refgenie_with_fasta.alias.resolve("t7") == digest


def _build_metadata_only_store(store_path: Path, fasta_path: Path, src_path: Path) -> str:
    """Create an on-disk store with collection metadata but NO sequence bytes.

    Mirrors what `genome init --store` produces via `_import_remote_collection` ->
    `local_store.add_sequence_collection(collection)`. Returns the collection digest.
    """
    from gtars.refget import RefgetStore

    # Source store built from a real FASTA: has both metadata and sequence bytes.
    src = RefgetStore.on_disk(str(src_path))
    meta, _ = src.add_sequence_collection_from_fasta(str(fasta_path))
    digest = meta.digest
    collection = src.get_collection(digest)

    # Target store gets metadata only -- no sequence bytes written.
    tgt = RefgetStore.on_disk(str(store_path))
    tgt.add_sequence_collection(collection)
    tgt.add_collection_alias("refgenie", "meta_only", digest)
    return digest


class TestStoreInitMetadataOnly:
    """A metadata-only store must not poison the store property (component tier).

    A `genome init --store <url>` imports collection metadata only (no sequence
    bytes). The store property used to eagerly load every sequence's bytes, raising
    OSError on such a store and poisoning every command that touched it. The fix
    removed eager loading; opening/resolving/listing must all succeed.
    """

    pytestmark = pytest.mark.component

    def test_metadata_only_store_does_not_poison(self, tmp_path, fixtures_path):
        genome_folder = tmp_path / "genomes"
        store_path = genome_folder / ".refget_store"
        store_path.parent.mkdir(parents=True, exist_ok=True)

        digest = _build_metadata_only_store(
            store_path, fixtures_path / "rCRSd.fa", tmp_path / "src"
        )

        r = Refgenie(database_engine=make_engine(), suppress_migrations=True)
        r.init(genome_folder=genome_folder)

        # The exact crash site: opening the store used to eagerly load bytes.
        store = r.refget_store
        assert store is not None

        # Alias resolution reads collection metadata -- never bytes.
        assert r.alias.resolve("meta_only") == digest

        # Listing genomes must not touch sequence bytes either.
        r.genome.list_all()


class TestApplyFhr:
    """apply_fhr upserts genome metadata columns and writes the store FHR sidecar.

    This is the single FHR -> catalog mapping both CLI callers (`genome init
    --fhr`, `genome set-metadata --fhr`) and the registry's post-build
    apply_metadata step funnel through, so its behavior is locked here.
    """

    pytestmark = pytest.mark.component

    def _init(self, refgenie, fixtures_path):
        return refgenie.genome.initialize_genome(
            fasta_file_path=fixtures_path / "rCRSd.fa",
            description="",
            alias_names=["myg"],
        )

    @pytest.mark.parametrize(
        "fhr, expected",
        [
            (
                {"genome": "Homo sapiens", "documentation": "A test genome.", "version": "myg"},
                {"species_name": "Homo sapiens", "description": "A test genome."},
            ),
            (
                {
                    "genome": "Homo sapiens",
                    "commonName": "human",
                    "taxon": {
                        "name": "Homo sapiens",
                        "uri": "https://identifiers.org/taxonomy:9606",
                    },
                    "assemblySource": "NCBI",
                    "accessionID": {"name": "GCA_000001405.15"},
                    "assemblyLevel": "chromosome",
                },
                {
                    "common_name": "human", "taxon_id": 9606, "assembly_source": "NCBI",
                    "assembly_accession": "GCA_000001405.15", "assembly_level": "chromosome",
                },
            ),
        ],
        ids=["basic-metadata", "taxon-and-assembly"],
    )
    def test_apply_fhr_populates_columns(self, refgenie_with_fasta, fixtures_path, fhr, expected):
        digest, _ = self._init(refgenie_with_fasta, fixtures_path)
        refgenie_with_fasta.genome.apply_fhr(digest, fhr)

        g = refgenie_with_fasta.genome.get(digest)
        for attr, value in expected.items():
            assert getattr(g, attr) == value, f"column {attr!r} not applied"
        # The record round-trips through the RefgetStore sidecar the server serves.
        assert refgenie_with_fasta.refget_store.get_fhr_metadata(digest).genome == fhr["genome"]

    @pytest.mark.parametrize(
        "first_fhr, expected",
        [
            (
                {"genome": "Homo sapiens", "documentation": "desc"},
                {"species_name": "Homo sapiens", "description": "desc"},
            ),
            (
                {"assemblySource": "UCSC", "accessionID": {"name": "GCA_1"}},
                {"assembly_source": "UCSC", "assembly_accession": "GCA_1"},
            ),
        ],
        ids=["species-and-description", "assembly-columns"],
    )
    def test_apply_fhr_partial_record_does_not_clobber(
        self, refgenie_with_fasta, fixtures_path, first_fhr, expected
    ):
        digest, _ = self._init(refgenie_with_fasta, fixtures_path)
        refgenie_with_fasta.genome.apply_fhr(digest, first_fhr)
        # A minimal record (a genome with no YAML) must leave existing columns intact.
        refgenie_with_fasta.genome.apply_fhr(digest, {"name": "myg"})
        g = refgenie_with_fasta.genome.get(digest)
        for attr, value in expected.items():
            assert getattr(g, attr) == value, f"column {attr!r} not applied"

    def test_missing_genome_raises(self, refgenie_with_fasta):
        with pytest.raises(MissingGenomeError):
            refgenie_with_fasta.genome.apply_fhr(
                "ffffffff1234567890abcdef12345678", {"genome": "X"}
            )


@pytest.mark.component
class TestGenomeGetMetadata:
    def test_returns_expected_keys_and_values(self, refgenie_with_genome):
        rgc = refgenie_with_genome
        digest = rgc.alias.resolve("rCRSd")
        metadata = rgc.genome.get_metadata(digest)
        assert metadata["digest"] == digest
        assert metadata["n_sequences"] >= 1
        assert metadata["total_length"] > 0
        assert metadata["source"] == "local"

    def test_unknown_digest_raises(self, refgenie_with_genome):
        from refgenie.exceptions import MissingGenomeError

        with pytest.raises(MissingGenomeError):
            refgenie_with_genome.genome.get_metadata("nonexistent_digest_00000000000000")


# ---------------------------------------------------------------------------
# Refgenie.init(): folder creation, idempotency, and empty config (unit tier)
#
# Covers the regression where REFGENIE_GENOME_STAGE_FOLDER was not auto-created
# at init.
# ---------------------------------------------------------------------------


def test_init_creates_genome_and_stage_folders(engine, tmp_path):
    genome_folder = tmp_path / "g"
    stage_folder = tmp_path / "s"
    assert not genome_folder.exists()
    assert not stage_folder.exists()

    r = Refgenie(database_engine=engine, suppress_migrations=True)
    r.init(genome_folder=genome_folder, genome_stage_folder=stage_folder)

    assert genome_folder.is_dir()
    assert stage_folder.is_dir()


def test_init_idempotent_when_folders_exist(engine, tmp_path):
    genome_folder = tmp_path / "g"
    stage_folder = tmp_path / "s"
    genome_folder.mkdir(parents=True)
    stage_folder.mkdir(parents=True)

    r = Refgenie(database_engine=engine, suppress_migrations=True)
    # Should not raise
    r.init(genome_folder=genome_folder, genome_stage_folder=stage_folder)

    assert genome_folder.is_dir()
    assert stage_folder.is_dir()


def test_init_no_stage_folder_when_none(engine, tmp_path, monkeypatch):
    # Force config.genome_stage_folder to None so init() truly receives no stage folder.
    from refgenie import config as config_module

    monkeypatch.setattr(config_module.config, "genome_stage_folder", None)

    genome_folder = tmp_path / "g"
    stage_folder = tmp_path / "s"
    assert not genome_folder.exists()

    r = Refgenie(database_engine=engine, suppress_migrations=True)
    r.init(genome_folder=genome_folder, genome_stage_folder=None)

    assert genome_folder.is_dir()
    # With no stage folder configured, init() must not create one.
    assert not stage_folder.exists()


def test_init_backend_idempotent(engine, tmp_path):
    """Regression: init_backend is idempotent.

    The persistent build catalog (nightly registry builds) re-runs init on every
    run. Configuration.version is unique, so a naive re-insert raised IntegrityError
    -- logged as an alarming ERROR that also masked genuine init failures.
    init_backend now skips the insert when a Configuration row already exists.
    """
    r = Refgenie(database_engine=engine, suppress_migrations=True)
    genome_folder = tmp_path / "genomes"
    genome_folder.mkdir()

    # First init creates the Configuration row.
    assert r.init_backend(genome_folder=genome_folder) is True
    # Re-init on an already-initialized backend is a no-op success, not an
    # IntegrityError, and does not insert a duplicate row.
    assert r.init_backend(genome_folder=genome_folder) is True

    with r._database_session as session:
        assert len(session.exec(select(Configuration)).all()) == 1


def test_initialization_produces_empty_config(engine, tmp_path):
    """After init, no asset classes or recipes are registered."""
    r = Refgenie(database_engine=engine, suppress_migrations=True)
    r.init(genome_folder=tmp_path / "genomes")
    assert len(r.recipe.list_all()) == 0
    assert len(r.asset_class.list_all()) == 0


def test_bare_init_never_writes_under_home(engine, tmp_path):
    """A bare init() must not resolve into the user's real ~/.refgenie.

    ``check_for_db_migrations`` can call a bare ``init()`` on its own, so the
    default has to be safe rather than fatal. The conftest redirect plus the
    autouse ``_isolated_default_genome_folder`` fixture are what make it so.
    """
    r = Refgenie(database_engine=engine, suppress_migrations=True)
    r.init()
    assert Path.home() / ".refgenie" not in r.genome_folder.parents
    assert r.genome_folder != Path.home() / ".refgenie" / "genomes"


class TestRefgenieInit:
    """Refgenie constructor input coercion."""

    def test_accepts_str_database_config_path(self, tmp_path):
        """Refgenie(database_config_path=<str>) coerces to Path without AttributeError."""
        config_path = tmp_path / "refgenie_db_config.yaml"
        assert Refgenie(database_config_path=str(config_path), suppress_migrations=True) is not None
        assert Refgenie(database_config_path=config_path, suppress_migrations=True) is not None


# ---------------------------------------------------------------------------
# Alias manager (refgenie.managers.alias) and the alias symlink tree
#
# Two surfaces live here:
#
# * **TestAliasCRUD** (unit) -- catalog-level alias operations: resolve,
#   exists, list/filter, get_for_genome, table, and the add/set/remove
#   lifecycle.
# * **TestAliasSymlinkTree** (component) -- the name-addressed view of the
#   data directory. For every alias a genome has,
#   ``alias/<name>/<group>/<asset>/`` mirrors the asset directory with the
#   genome digest in each filename rewritten to that alias. Two load-bearing,
#   previously-wrong properties are pinned here: the tree is rooted at
#   ``Asset.path`` (not a seek key's parent), and each alias names its own
#   tree.
# ---------------------------------------------------------------------------


class TestAliasCRUD:
    """Catalog-level alias operations (unit tier)."""

    def test_resolve(self, refgenie_session):
        """resolve returns the 32-char seqcol digest by positional or keyword arg."""
        by_pos = refgenie_session.alias.resolve("rCRSd")
        by_kw = refgenie_session.alias.resolve(name="rCRSd")
        assert by_pos == by_kw
        assert len(by_pos) == 32
        with pytest.raises(MissingAliasError):
            refgenie_session.alias.resolve("nonexistent_alias_12345")

    def test_exists(self, refgenie_session):
        assert refgenie_session.alias.exists("rCRSd")
        assert not refgenie_session.alias.exists("nonexistent_alias_12345")

    def test_list_all_and_genome_filter(self, refgenie_session):
        """list_all contains rCRSd; filtered by digest is non-empty and all rows match."""
        digest = refgenie_session.alias.resolve("rCRSd")
        names = {a.name for a in refgenie_session.alias.list_all()}
        assert "rCRSd" in names

        filtered = list(refgenie_session.alias.list_all(genome_digest=digest))
        assert filtered
        assert all(a.genome_digest == digest for a in filtered)

    def test_get_for_genome(self, refgenie_session):
        """get_for_genome returns the alias names of a genome, including the queried one."""
        digest = refgenie_session.alias.resolve("rCRSd")
        assert "rCRSd" in refgenie_session.alias.get_for_genome(digest)

    def test_remove_missing_raises(self, refgenie_minimal):
        with pytest.raises(MissingAliasError):
            refgenie_minimal.alias.remove("nonexistent_alias_12345")

    def test_add_genome_with_alias(self, refgenie_minimal):
        """genome.add with an alias list registers the alias, resolvable back to the digest."""
        digest = "abcdef1234567890abcdef1234567890"
        genome = refgenie_minimal.genome.add(
            digest, "Test genome created with alias", ["test_genome_with_alias"]
        )
        assert genome.digest == digest
        assert refgenie_minimal.alias.exists("test_genome_with_alias")
        assert refgenie_minimal.alias.resolve("test_genome_with_alias") == digest

    def test_set_and_remove_alias(self, refgenie_minimal):
        """set_genome_alias then alias.remove leaves the alias gone."""
        digest = "bbbbbbbb1234567890abcdef12345678"
        refgenie_minimal.genome.add(digest, "Test genome for alias removal", [])
        refgenie_minimal.set_genome_alias("test_remove_genome_alias", digest)
        assert refgenie_minimal.alias.exists("test_remove_genome_alias")

        refgenie_minimal.alias.remove("test_remove_genome_alias")
        assert not refgenie_minimal.alias.exists("test_remove_genome_alias")

    def test_set_genome_alias_new_genome_refreshes_alias_tree(self, refgenie_minimal, monkeypatch):
        """Setting an alias for a genome the catalog has never seen must behave
        like every other exit path: register the alias exactly once and refresh
        the symlink tree exactly once (not register twice, not skip the tree)."""
        symlink_calls = []
        add_calls = []

        real_add = refgenie_minimal._alias_manager.add

        monkeypatch.setattr(
            refgenie_minimal,
            "_symlink_alias",
            lambda **kwargs: symlink_calls.append(kwargs) or [],
        )
        monkeypatch.setattr(
            refgenie_minimal._alias_manager,
            "add",
            lambda *args, **kwargs: (add_calls.append(args), real_add(*args, **kwargs))[1],
        )

        refgenie_minimal.set_genome_alias("brand_new", genome_digest="0" * 48)

        assert len(add_calls) == 1
        assert len(symlink_calls) == 1


# --- Symlink-tree scaffolding (component tier) ------------------------------

SECOND_ALIAS = "rCRSd_alt"
ALIAS_DIR = "alias"


def _alias_dir(r: Refgenie, alias: str, group: str, asset: str) -> Path:
    return r.genome_folder / ALIAS_DIR / alias / group / asset


class TestAliasSymlinkTree:
    """The name-addressed symlink view of the data directory (component tier)."""

    pytestmark = pytest.mark.component

    def test_get_asset_dir_is_the_recorded_asset_path(self, refgenie_fs, fixtures_path):
        """The asset directory is Asset.path resolved against the genome folder."""
        r = refgenie_fs
        asset = add_asset_from_files(
            r,
            asset_class_name="fasta",
            asset_group_name="fasta",
            asset_name="test",
            rel_dir="data/fasta_payload",
            files=["{genome}.fa", "{genome}.fa.fai", "{genome}.chrom.sizes"],
        )
        assert r.asset.get_asset_dir(
            genome_name=GENOME,
            asset_group_name="fasta",
            asset_name="test",
        ) == (r.genome_folder / asset.path)

    def test_get_asset_dir_rejects_an_incomplete_asset(self, refgenie_fs, monkeypatch):
        """An asset with no path has no directory; that is an error, not a guess."""
        r = refgenie_fs
        asset = add_asset_from_files(
            r,
            asset_class_name="fasta",
            asset_group_name="fasta",
            asset_name="test",
            rel_dir="data/fasta_payload",
            files=["{genome}.fa", "{genome}.fa.fai", "{genome}.chrom.sizes"],
        )
        # An incomplete asset cannot be inserted through the normal path, so stand
        # one in: what matters is that get_asset_dir refuses to build a path from a
        # null column rather than producing genome_folder itself.
        asset.path = None
        monkeypatch.setattr(r.asset, "get", lambda **kwargs: asset)

        with pytest.raises(ValueError, match="path is not set"):
            r.asset.get_asset_dir(
                genome_name=GENOME,
                asset_group_name="fasta",
                asset_name="test",
            )

    def test_nested_seek_key_value_does_not_truncate_the_tree(self, refgenie_fs, fixtures_path):
        """A seek key nested in a subdirectory must not become the symlink root.

        Regression test. Deriving the source from ``seek().parent`` rooted the tree
        at ``sub/``, so everything above it was silently omitted.
        """
        r = refgenie_fs
        r.asset_class.add(fixtures_path / "nested_seek_key_asset_class.yaml")
        add_asset_from_files(
            r,
            asset_class_name="nested_seek_key",
            asset_group_name="nested_seek_key",
            asset_name="test",
            rel_dir="data/nested_payload",
            files=["sub/inner/{genome}.idx", "{genome}.toplevel"],
        )

        alias_dir = _alias_dir(r, GENOME, "nested_seek_key", "test")
        # The nested seek key is still linked, at its real depth on disk.
        assert (alias_dir / "sub" / "inner" / f"{GENOME}.idx").is_symlink()
        # And the file that lives above it — invisible to the old .parent trick,
        # which rooted the tree at sub/ and never saw this file at all.
        assert (alias_dir / f"{GENOME}.toplevel").is_symlink()
        assert (alias_dir / f"{GENOME}.toplevel").resolve().read_text() == "payload\n"

    def test_non_path_default_seek_key_still_gets_a_tree(self, refgenie_fs, fixtures_path):
        """An asset with no path seek keys is still linkable.

        The old code raised TypeError here, denying a name-addressed view because
        of an unrelated property of one of the asset's seek keys.
        """
        r = refgenie_fs
        r.asset_class.add(fixtures_path / "metadata_only_asset_class.yaml")
        add_asset_from_files(
            r,
            asset_class_name="metadata_only",
            asset_group_name="metadata_only",
            asset_name="test",
            rel_dir="data/metadata_payload",
            files=["{genome}.notes"],
            custom_seek_keys={"software_version": "1.2.3", "build_info": "{}"},
        )

        alias_dir = _alias_dir(r, GENOME, "metadata_only", "test")
        assert (alias_dir / f"{GENOME}.notes").is_symlink()

    def test_each_alias_tree_is_named_for_its_own_alias(self, refgenie_fs, fixtures_path):
        """Regression test: a second alias must not inherit the first alias's names."""
        r = refgenie_fs
        digest = genome_digest(r)
        r.alias.add(name=SECOND_ALIAS, genome_digest=digest)
        add_asset_from_files(
            r,
            asset_class_name="fasta",
            asset_group_name="fasta",
            asset_name="test",
            rel_dir="data/fasta_payload",
            files=["{genome}.fa", "{genome}.fa.fai", "{genome}.chrom.sizes"],
        )
        r._symlink_alias(alias_name=GENOME, asset_group_name="fasta", asset_name="test")

        for alias in (GENOME, SECOND_ALIAS):
            alias_dir = _alias_dir(r, alias, "fasta", "test")
            assert (alias_dir / f"{alias}.fa").is_symlink(), f"missing payload for {alias}"
            assert not (alias_dir / f"{digest}.fa").exists(), f"digest survived in {alias} tree"

        # The decisive assertion: neither tree is named for the other's alias.
        assert not (_alias_dir(r, SECOND_ALIAS, "fasta", "test") / f"{GENOME}.fa").exists()

    def test_add_with_genome_digest_still_builds_the_alias_tree(self, refgenie_fs):
        """``genome_digest=`` is documented as equivalent to ``genome_name=``,
        so add() must build the alias symlink tree either way."""
        r = refgenie_fs
        digest = genome_digest(r)
        rel_dir = "data/fasta_payload"
        asset_dir = r.genome_folder / rel_dir
        asset_dir.mkdir(parents=True, exist_ok=True)
        for name in (f"{digest}.fa", f"{digest}.fa.fai", f"{digest}.chrom.sizes"):
            (asset_dir / name).write_text("payload\n")

        r.add(
            genome_digest=digest,
            asset_group_name="fasta",
            asset_name="test",
            path=Path(rel_dir),
            asset_class_name="fasta",
        )

        alias_dir = _alias_dir(r, GENOME, "fasta", "test")
        assert alias_dir.is_dir(), f"no alias tree at {alias_dir}"
        assert (alias_dir / f"{GENOME}.fa").is_symlink()

    def test_group_level_symlink_refresh_without_an_asset_name(self, refgenie_fs):
        """``_symlink_alias(group)`` with no asset name must resolve the group's
        default asset (by name, not digest) and refresh its tree."""
        r = refgenie_fs
        add_asset_from_files(
            r,
            asset_class_name="fasta",
            asset_group_name="fasta",
            asset_name="test",
            rel_dir="data/fasta_payload",
            files=["{genome}.fa", "{genome}.fa.fai", "{genome}.chrom.sizes"],
        )

        r._symlink_alias(alias_name=GENOME, asset_group_name="fasta")

        assert _alias_dir(r, GENOME, "fasta", "test").is_dir()

    def test_repeated_symlinking_is_idempotent_and_repoints(self, tmp_path):
        """A second run replaces stale links rather than silently keeping them."""
        genome_digest = "abc123"
        first_src = tmp_path / "first"
        first_src.mkdir()
        (first_src / f"{genome_digest}.fa").write_text("first\n")

        target = tmp_path / ALIAS_DIR / "hg38"
        mapping = {"hg38": target}

        create_alias_symlinks(
            src_path=first_src, target_paths_mapping=mapping, genome_digest=genome_digest
        )
        link = target / "hg38.fa"
        assert link.resolve().read_text() == "first\n"

        second_src = tmp_path / "second"
        second_src.mkdir()
        (second_src / f"{genome_digest}.fa").write_text("second\n")

        create_alias_symlinks(
            src_path=second_src, target_paths_mapping=mapping, genome_digest=genome_digest
        )
        assert link.is_symlink()
        assert link.resolve().read_text() == "second\n"

    def test_a_real_file_in_the_alias_tree_is_not_destroyed(self, tmp_path):
        """A non-symlink at a destination is data, not a link we own. Raise."""
        genome_digest = "abc123"
        src = tmp_path / "src"
        src.mkdir()
        (src / f"{genome_digest}.fa").write_text("payload\n")

        target = tmp_path / ALIAS_DIR / "hg38"
        target.mkdir(parents=True)
        squatter = target / "hg38.fa"
        squatter.write_text("real data\n")

        with pytest.raises(FileExistsError, match="Refusing to replace non-symlink"):
            create_alias_symlinks(
                src_path=src,
                target_paths_mapping={"hg38": target},
                genome_digest=genome_digest,
            )
        assert squatter.read_text() == "real data\n"


# ---------------------------------------------------------------------------
# getseq() using the RefgetStore backend (unit tier)
#
# Tests use in-memory RefgetStore to avoid needing samtools or FASTA files on
# disk.
# ---------------------------------------------------------------------------

RCRSD_LENGTH = 33138


def _local_refgenie_with_bytes(tmp_path: Path, fasta_path: Path):
    """A LocalMode Refgenie whose on-disk store has full sequence bytes.

    The store is written to disk *before* Refgenie opens it, so Refgenie opens
    it lazily (no eager sequence load) -- exercising the on-demand load path.

    Returns (refgenie, collection_digest, sequence_name).
    """
    from gtars.refget import RefgetStore

    genome_folder = tmp_path / "genomes"
    store_path = genome_folder / ".refget_store"
    store_path.parent.mkdir(parents=True, exist_ok=True)

    store = RefgetStore.on_disk(str(store_path))
    meta, _ = store.add_sequence_collection_from_fasta(str(fasta_path))
    digest = meta.digest
    store.add_collection_alias("refgenie", "rCRSd", digest)

    r = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    r.init(genome_folder=genome_folder)
    return r, digest, "rCRSd"


def test_getseq_whole_chromosome_lazy_store(tmp_path, fixtures_path):
    """Whole-chromosome getseq loads bytes on demand from a lazily-opened store.

    With eager load_all_sequences() gone, the record returned by
    get_sequence_by_name has no bytes until explicitly loaded, so record.decode()
    initially returns None. getseq must load the sequence and re-fetch.
    """
    r, _digest, name = _local_refgenie_with_bytes(tmp_path, fixtures_path / "rCRSd.fa")

    seq = r.getseq("rCRSd", name)
    assert isinstance(seq, str)
    assert len(seq) == RCRSD_LENGTH


def test_getseq_substring_lazy_store(tmp_path, fixtures_path):
    """Substring getseq works against a lazily-opened store (explicit load path)."""
    r, _digest, name = _local_refgenie_with_bytes(tmp_path, fixtures_path / "rCRSd.fa")

    whole = r.getseq("rCRSd", name)
    sub = r.getseq("rCRSd", f"{name}:0-10")
    assert len(sub) == 10
    assert whole.startswith(sub)


def test_getseq_metadata_only_routes_to_remote(tmp_path, fixtures_path):
    """getseq on a metadata-only (store-initialized) genome routes to remote fallback.

    A metadata-only collection is present in the store, so get_sequence_by_name
    returns a bytes-less record (NOT a KeyError). The fix must recognize the
    missing local bytes and fall back to the genome's remote_url rather than
    raising OSError -- this is the getseq facet of the store-init poisoning bug.
    """
    from gtars.refget import RefgetStore

    from refgenie.db.tables import Genome, Alias

    genome_folder = tmp_path / "genomes"
    store_path = genome_folder / ".refget_store"
    store_path.parent.mkdir(parents=True, exist_ok=True)

    # Source store with real bytes, used only to produce a valid collection.
    src = RefgetStore.on_disk(str(tmp_path / "src"))
    meta, _ = src.add_sequence_collection_from_fasta(str(fixtures_path / "rCRSd.fa"))
    digest = meta.digest

    # Local store: metadata only (no sequence bytes) -- the store-init state.
    tgt = RefgetStore.on_disk(str(store_path))
    tgt.add_sequence_collection(src.get_collection(digest))
    tgt.add_collection_alias("refgenie", "meta_only", digest)

    r = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    r.init(genome_folder=genome_folder)

    # Genome row with a remote_url so the fallback has a URL (as store-init records).
    with r._database_session as session:
        session.add(Genome(digest=digest, description="meta", remote_url="https://x"))
        session.add(Alias(name="meta_only", genome_digest=digest))
        session.commit()

    mock_record = MagicMock()
    mock_record.decode.return_value = "ACGTACGT"
    mock_record.metadata.sha512t24u = "fake_digest"
    mock_record.metadata.length = 8

    with patch.object(r, "_fetch_remote_sequence", return_value=mock_record) as mock_fetch:
        result = r.getseq("meta_only", "rCRSd")

    # The bytes-less local record must NOT raise OSError; it routes to remote.
    mock_fetch.assert_called_once()
    assert result == "ACGTACGT"


class TestParseLocus:
    """Tests for the _parse_locus helper function."""

    @pytest.mark.parametrize(
        "locus, expected",
        [
            ("chr1", ("chr1", None, None)),  # name only
            ("chr1:0-1000", ("chr1", 0, 1000)),  # name with range
            ("chr1:500", ("chr1", 500, None)),  # name with start only
            ("V01146.1:0-10", ("V01146.1", 0, 10)),  # dot in name
            ("ENST00000-1:0-10", ("ENST00000-1", 0, 10)),  # hyphen in name
            ("gi|12345|ref|NC_001:0-100", ("gi|12345|ref|NC_001", 0, 100)),  # pipes in name
            ("invalid locus", ValueError),  # space is invalid
            ("", ValueError),  # empty is invalid
        ],
    )
    def test_parse_locus(self, locus, expected):
        """_parse_locus splits name/start/end across name char-classes; rejects bad input."""
        if expected is ValueError:
            with pytest.raises(ValueError, match="Invalid locus format"):
                _parse_locus(locus)
        else:
            assert _parse_locus(locus) == expected
