"""
Tests for the genome manager (refgenie.managers.genome / sources.genomes) and the
genome init/build flow.

Unit-tier classes cover catalog queries, removal cascade, build-input validation (pure static-method checks), and the genome-init
failure/exit contract (via lightweight doubles). Component-tier classes exercise
the real initialize_and_build flow, alias self-repair on re-init, and the
metadata-only store regression on a real on-disk gtars store.

Also here: the alias manager and its symlink tree. Database init lives in
test_database.py and sequence retrieval in test_sequence.py; the ORM get
round-trip stays in test_db.py.
"""

from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from refgenie import Refgenie
from refgenie.cli.commands.genome import GenomeInitModel, handle_genome_init
from refgenie.exceptions import MissingAliasError, MissingGenomeError
from refgenie.utils.symlinks import create_alias_symlinks
from refgenie.models import GenomeAlias, GenomeDigest
from tests.helpers import GENOME, add_asset_from_files, fake_digest, genome_digest, make_engine


class TestGenomeIdentifierTypes:
    """GenomeDigest and GenomeAlias: str subclasses that say which kind of
    genome identifier a value is (unit tier)."""

    def test_digest_accepts_a_real_digest(self, refgenie_session):
        real = refgenie_session.alias.resolve("rCRSd")
        digest = GenomeDigest(str(real))
        assert digest == real
        assert isinstance(digest, str) and isinstance(digest, GenomeDigest)

    @pytest.mark.parametrize(
        "value", ["", "hg38", "a" * 31, "a" * 33, "a" * 31 + "/", "a" * 31 + "="]
    )
    def test_digest_rejects_a_malformed_value(self, value):
        with pytest.raises(ValueError, match="not a genome digest"):
            GenomeDigest(value)

    @pytest.mark.parametrize("value", ["", "hg38/fasta"])
    def test_alias_rejects_empty_or_slash(self, value):
        with pytest.raises(ValueError, match="not a genome alias"):
            GenomeAlias(value)

    def test_alias_allows_names_and_qualified_names(self):
        assert GenomeAlias("hg38") == "hg38"
        assert GenomeAlias("ucsc::hg38") == "ucsc::hg38"
        assert isinstance(GenomeAlias("hg38"), str)

    def test_pydantic_fields_validate_and_produce_the_subclass(self):
        from pydantic import BaseModel, ValidationError

        class Ref(BaseModel):
            digest: GenomeDigest
            alias: GenomeAlias

        ref = Ref(digest="a" * 32, alias="hg38")
        assert type(ref.digest) is GenomeDigest
        assert type(ref.alias) is GenomeAlias
        assert ref.model_dump() == {"digest": "a" * 32, "alias": "hg38"}
        with pytest.raises(ValidationError):
            Ref(digest="hg38", alias="hg38")
        with pytest.raises(ValidationError):
            Ref(digest="a" * 32, alias="")

    def test_fake_digest_is_deterministic_and_well_formed(self):
        assert fake_digest("x") == fake_digest("x") != fake_digest("y")
        assert isinstance(fake_digest("x"), GenomeDigest)


class TestGenomeCatalog:
    """Catalog-level genome queries (unit tier)."""

    def test_list_all_and_exists(self, refgenie_session):
        """The built genome is listed and exists() agrees; unknown digests do not."""
        digest = refgenie_session.alias.resolve("rCRSd")
        assert digest in {g.digest for g in refgenie_session.genome.list_all()}
        assert refgenie_session.genome.exists(digest)
        assert not refgenie_session.genome.exists("nonexistent_genome")

    def test_alias_resolve_returns_a_genome_digest(self, refgenie_session):
        """alias.resolve maps an alias to its digest, typed as a GenomeDigest."""
        digest = refgenie_session.alias.resolve(GenomeAlias("rCRSd"))
        assert isinstance(digest, GenomeDigest)
        assert digest != "rCRSd"

    @pytest.mark.parametrize("backend", ["federated", "store", "sql"])
    def test_alias_resolve_has_no_digest_fallback(self, refgenie_minimal, backend):
        """A digest is not an alias: every backend refuses to resolve one."""
        from refgenie.managers.alias import AliasManager

        r = refgenie_minimal
        digest = fake_digest("no-fallback")
        r.genome.add(digest, "a genome", [])
        aliases = {
            "federated": r.alias,
            "store": r.alias._local,
            "sql": AliasManager(r.database_engine),
        }[backend]
        aliases.add(GenomeAlias("named"), digest)
        resolved = aliases.resolve(GenomeAlias("named"))
        assert resolved == digest and isinstance(resolved, GenomeDigest)
        with pytest.raises(MissingAliasError):
            aliases.resolve(digest)

    def test_alias_spelled_like_a_digest_is_an_alias(self, refgenie_minimal):
        """An alias spelled like another genome's digest names its own genome.

        Aliases resolve only as aliases, so the shadowing case cannot be
        confused: the digest-shaped alias leads to the genome it was set on.
        """
        r = refgenie_minimal
        genome_a, genome_b = fake_digest("a"), fake_digest("b")
        r.genome.add(genome_a, "a", [])
        r.genome.add(genome_b, "b", [])
        r.alias.add(GenomeAlias(genome_b), genome_a)
        assert r.alias.resolve(GenomeAlias(genome_b)) == genome_a
        assert r.genome.get(genome_b).digest == genome_b


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
# Regression: `build.run` signals failure by RETURNING None rather than
# raising, so `handle_genome_init` wrapping the call in
# try/except caught nothing -- a failed build logged "built successfully" and
# exited 0. The orchestrator gates on that exit status, so it recorded the genome
# as initialized and the real error surfaced later as a misleading "not found".


def _init_cmd(**overrides):
    """Stand-in for the parsed `genome init --fasta ... --name ...` command."""
    fields = dict(
        name=["athaliana"],
        fasta="/tmp/does-not-matter.fa",
        server=None,
        store=None,
        namespace=None,
        digest=None,
        description="",
        species=None,
        fhr=None,
        force=False,
        build=True,
    )
    fields.update(overrides)
    return SimpleNamespace(**fields)


# The double must not drift from the real parsed model.
assert set(_init_cmd().__dict__) <= set(GenomeInitModel.model_fields)


def _init_refgenie(build_result):
    """A refgenie double whose build.run outcome the caller controls."""
    rg = MagicMock()
    rg.genome.initialize_genome.return_value = ("AHyVecap0miRdrz6ZffizkZQUOJJBO2c", True)
    rg.recipe.exists.return_value = True
    if isinstance(build_result, Exception):
        rg.build.run.side_effect = build_result
    else:
        rg.build.run.return_value = build_result
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
        refgenie.build.run.assert_called_once()
        assert "built successfully" not in caplog.text

    def test_succeeds_when_build_returns_asset(self):
        """A genuine success must NOT exit; that is what makes the failure cases meaningful."""
        refgenie = _init_refgenie(object())
        handle_genome_init(_init_cmd(), refgenie)
        refgenie.build.run.assert_called_once()


@pytest.fixture
def init_t7(refgenie_with_fasta, fixtures_path):
    """Run build.initialize_and_build for the 't7' genome from the rCRSd fixture FASTA."""

    def _init(**overrides):
        kwargs = dict(
            fasta_file_path=fixtures_path / "rCRSd.fa",
            genome_names=["t7"],
            description="Test genome",
        )
        kwargs.update(overrides)
        return refgenie_with_fasta.build.initialize_and_build(**kwargs)

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
    """The build.initialize_and_build convenience method (component tier)."""

    pytestmark = pytest.mark.component

    def test_init_with_fasta_auto_builds(self, refgenie_with_fasta, init_t7):
        """build.initialize_and_build creates both the genome and the fasta asset."""
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
        refgenie_with_fasta.build.run(
            recipe_name="fasta",
            genome_alias="t7",
            asset_group_name="fasta",
        )
        asset = refgenie_with_fasta.asset.get(
            genome_digest=genome_digest,
            asset_group_name="fasta",
            asset_name="default",
        )
        assert asset is not None

    def test_init_idempotency_with_existing(self, init_t7):
        """Calling build.initialize_and_build twice with use_existing=True does not crash."""
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
    bytes). The store property must not eagerly load every sequence's bytes:
    otherwise it raises OSError on such a store and poisons every command that
    touches it. Opening, resolving and listing must all succeed.
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
        r.database.init(genome_folder=genome_folder)

        # Opening the store must not load sequence bytes.
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
                    "common_name": "human",
                    "taxon_id": 9606,
                    "assembly_source": "NCBI",
                    "assembly_accession": "GCA_000001405.15",
                    "assembly_level": "chromosome",
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
#   genome digest in each filename rewritten to that alias. Two load-bearing
#   properties are pinned here: the tree is rooted at
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

    def test_get_for_genomes_matches_per_genome_lookup(self, refgenie_session):
        """get_for_genomes agrees with get_for_genome for every digest, and maps
        a digest with no aliases to an empty list instead of raising."""
        alias = refgenie_session.alias
        digests = sorted({a.genome_digest for a in alias.list_all()})
        unknown = GenomeDigest("0" * 32)

        batched = alias.get_for_genomes([*digests, unknown])

        for digest in digests:
            assert sorted(batched[digest]) == sorted(alias.get_for_genome(digest))
        assert batched[unknown] == []
        assert alias.get_for_genomes([]) == {}

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

        real_add = refgenie_minimal.alias.add

        monkeypatch.setattr(
            refgenie_minimal.asset.tree,
            "render_genome",
            lambda *args, **kwargs: symlink_calls.append((args, kwargs)),
        )
        monkeypatch.setattr(
            refgenie_minimal.alias,
            "add",
            lambda *args, **kwargs: (add_calls.append(args), real_add(*args, **kwargs))[1],
        )

        refgenie_minimal.set_genome_alias("brand_new", genome_digest="0" * 48)

        assert len(add_calls) == 1
        assert len(symlink_calls) == 1


# --- Symlink-tree scaffolding (component tier) ------------------------------

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
            genome_digest=genome_digest(r),
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
                genome_digest=genome_digest(r),
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
        # And the file that lives above it — a tree rooted at a seek key's
        # parent (sub/) would never see this file at all.
        assert (alias_dir / f"{GENOME}.toplevel").is_symlink()
        assert (alias_dir / f"{GENOME}.toplevel").resolve().read_text() == "payload\n"

    def test_non_path_default_seek_key_still_gets_a_tree(self, refgenie_fs, fixtures_path):
        """An asset with no path seek keys is still linkable.

        A seek key that is not a path must not raise TypeError here: that would
        deny a name-addressed view because of an unrelated property of one of
        the asset's seek keys.
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

    def test_add_with_a_digest_builds_the_alias_tree(self, refgenie_fs):
        """add() takes a digest and builds the alias symlink tree for every
        alias of the genome."""
        r = refgenie_fs
        digest = genome_digest(r)
        rel_dir = "data/fasta_payload"
        asset_dir = r.genome_folder / rel_dir
        asset_dir.mkdir(parents=True, exist_ok=True)
        for name in (f"{digest}.fa", f"{digest}.fa.fai", f"{digest}.chrom.sizes"):
            (asset_dir / name).write_text("payload\n")

        r.asset.content.add(
            genome_digest=digest,
            asset_group_name="fasta",
            asset_name="test",
            path=Path(rel_dir),
            asset_class_name="fasta",
        )

        alias_dir = _alias_dir(r, GENOME, "fasta", "test")
        assert alias_dir.is_dir(), f"no alias tree at {alias_dir}"
        assert (alias_dir / f"{GENOME}.fa").is_symlink()

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
