"""Federated multi-store serving: store registry, alias policy, and sync.

These build real on-disk RefgetStores and a real catalog, so they are
``component`` tier.
"""

import logging

import pytest

from refgenie import Refgenie
from refgenie.db.tables import StoreType
from refgenie.exceptions import (
    FederatedAliasError,
    MissingAliasError,
    MissingStoreError,
    StoreExistsError,
)
from refgenie.managers.alias import AliasManager

from .helpers import make_engine

pytestmark = pytest.mark.component


def _server_rgc():
    rg = Refgenie(database_engine=make_engine(), suppress_migrations=True, server_mode=True)
    rg._create_db_and_tables()
    return rg


def _build_on_disk_store(path, fasta, alias, namespace="name"):
    """Build an on-disk RefgetStore from a FASTA with one curated alias."""
    from refget.store import RefgetStore

    store = RefgetStore.on_disk(str(path))
    meta, _ = store.add_sequence_collection_from_fasta(str(fasta))
    store.add_collection_alias(namespace, alias, meta.digest)
    return meta.digest


# ---------------------------------------------------------------------------
# StoreManager CRUD
# ---------------------------------------------------------------------------


class TestStoreRegistry:
    def test_add_list_remove(self):
        rg = _server_rgc()
        rg.store.add("jungle", "url://j", StoreType.remote, priority=10)
        rg.store.add("vgp", "url://v", StoreType.remote, priority=20)
        assert [s.name for s in rg.store.enabled_stores()] == ["jungle", "vgp"]
        rg.store.remove("vgp")
        assert [s.name for s in rg.store.enabled_stores()] == ["jungle"]

    def test_add_duplicate_raises(self):
        rg = _server_rgc()
        rg.store.add("jungle", "url://j", StoreType.remote)
        with pytest.raises(StoreExistsError):
            rg.store.add("jungle", "url://other", StoreType.remote)

    def test_remove_missing_raises(self):
        rg = _server_rgc()
        with pytest.raises(MissingStoreError):
            rg.store.remove("ghost")

    def test_enabled_stores_priority_order(self):
        rg = _server_rgc()
        rg.store.add("c", "url://c", StoreType.remote, priority=30)
        rg.store.add("a", "url://a", StoreType.remote, priority=10)
        rg.store.add("b", "url://b", StoreType.remote, priority=20)
        assert [s.name for s in rg.store.enabled_stores()] == ["a", "b", "c"]


# ---------------------------------------------------------------------------
# Alias collision policy
# ---------------------------------------------------------------------------


class TestAliasPriority:
    def _setup(self, rg):
        """Two stores (jungle prio 10, vgp prio 20) with two genomes."""
        rg.store.add("jungle", "url://j", StoreType.remote, priority=10)
        rg.store.add("vgp", "url://v", StoreType.remote, priority=20)
        rg.genome.add("A" * 32, "genome A", [], store_name="jungle")
        rg.genome.add("B" * 32, "genome B", [], store_name="vgp")
        return AliasManager(rg.database_engine), {"jungle": 10, "vgp": 20}

    def test_higher_priority_wins_bare_name(self, caplog):
        rg = _server_rgc()
        alias, priorities = self._setup(rg)

        # vgp registers 'hg38' -> B first; then jungle registers 'hg38' -> A.
        alias.federated_sync("hg38", "B" * 32, "vgp", priorities)
        with caplog.at_level(logging.WARNING):
            outcome = alias.federated_sync("hg38", "A" * 32, "jungle", priorities)

        assert outcome == "repointed"
        # Bare name now resolves to the higher-priority store's digest.
        assert alias.resolve("hg38") == "A" * 32
        # The loser is reachable by its qualified name.
        assert alias.resolve("vgp::hg38") == "B" * 32
        # The winner is reachable by its qualified name too.
        assert alias.resolve("jungle::hg38") == "A" * 32
        # A WARNING named both stores.
        assert any("collision on 'hg38'" in r.message for r in caplog.records)

    def test_lower_priority_loses_bare_name(self, caplog):
        rg = _server_rgc()
        alias, priorities = self._setup(rg)

        # jungle registers first (wins), then vgp offers the same name.
        alias.federated_sync("hg38", "A" * 32, "jungle", priorities)
        with caplog.at_level(logging.WARNING):
            outcome = alias.federated_sync("hg38", "B" * 32, "vgp", priorities)

        assert outcome == "kept"
        assert alias.resolve("hg38") == "A" * 32
        assert alias.resolve("vgp::hg38") == "B" * 32
        assert any("collision on 'hg38'" in r.message for r in caplog.records)

    def test_same_digest_dedups(self):
        rg = _server_rgc()
        alias, priorities = self._setup(rg)
        assert alias.federated_sync("shared", "A" * 32, "jungle", priorities) == "inserted"
        # Same (alias, digest) from another store is a natural no-op.
        rg.genome.add("A" * 32, "genome A", [], store_name="vgp") if False else None
        assert alias.federated_sync("shared", "A" * 32, "vgp", priorities) == "noop"
        assert alias.resolve("shared") == "A" * 32

    def test_qualified_unknown_raises(self):
        rg = _server_rgc()
        alias, priorities = self._setup(rg)
        alias.federated_sync("hg38", "A" * 32, "jungle", priorities)
        with pytest.raises(MissingAliasError):
            alias.resolve("nostore::hg38")


# ---------------------------------------------------------------------------
# End-to-end: register + sync two on-disk stores, resolve across them
# ---------------------------------------------------------------------------


class TestStoreSyncIntegration:
    def test_two_stores_sync_dedup_and_route(self, tmp_path, fixtures_path):
        from refgenie.cli.commands.store import _sync_one_store

        # Two stores share the alias 'shared' for DIFFERENT genomes, and each
        # has one unique genome.
        jungle_shared = _build_on_disk_store(
            tmp_path / "jungle", fixtures_path / "rCRSd.fa", "shared"
        )
        vgp_shared = _build_on_disk_store(tmp_path / "vgp", fixtures_path / "t7.fa", "shared")

        rg = _server_rgc()
        rg.store.add("jungle", str(tmp_path / "jungle"), StoreType.on_disk, priority=10)
        rg.store.add("vgp", str(tmp_path / "vgp"), StoreType.on_disk, priority=20)

        priorities = {s.name: s.priority for s in rg.store.list_all()}
        for srow in rg.store.enabled_stores():
            _registered, failures = _sync_one_store(rg, srow, priorities, 100)
            assert failures == 0

        # Both genomes registered, each owned by its store.
        assert rg.genome.get(jungle_shared).store_name == "jungle"
        assert rg.genome.get(vgp_shared).store_name == "vgp"

        # Alias tie resolved by priority; loser reachable by qualified name.
        assert rg.alias.resolve("shared") == jungle_shared
        assert rg.alias.resolve("vgp::shared") == vgp_shared

        # The router routes each genome to its owning store.
        jungle_store = rg.store_router.store_for_genome(jungle_shared, store_name="jungle")
        assert jungle_store.is_collection_loaded(jungle_shared)
        vgp_store = rg.store_router.store_for_genome(vgp_shared, store_name="vgp")
        assert vgp_store.is_collection_loaded(vgp_shared)

    def test_identical_genome_across_stores_dedups(self, tmp_path, fixtures_path):
        """The same FASTA in two stores yields one Genome row (digest dedup)."""
        from refgenie.cli.commands.store import _sync_one_store

        digest_a = _build_on_disk_store(tmp_path / "a", fixtures_path / "rCRSd.fa", "rCRSd")
        digest_b = _build_on_disk_store(tmp_path / "b", fixtures_path / "rCRSd.fa", "rCRSd")
        assert digest_a == digest_b  # content-addressed: identical

        rg = _server_rgc()
        rg.store.add("a", str(tmp_path / "a"), StoreType.on_disk, priority=10)
        rg.store.add("b", str(tmp_path / "b"), StoreType.on_disk, priority=20)
        priorities = {s.name: s.priority for s in rg.store.list_all()}
        for srow in rg.store.enabled_stores():
            _sync_one_store(rg, srow, priorities, 100)

        # One genome row, owned by the higher-priority store.
        genomes = list(rg.genome.list_all())
        assert len(genomes) == 1
        assert genomes[0].store_name == "a"

    def test_resync_backfills_fhr_columns(self, tmp_path, fixtures_path):
        """A sidecar added to the store after the first sync lands on re-sync."""
        from refget.store import RefgetStore

        from refgenie.cli.commands.store import _sync_one_store
        from refgenie.managers.genome import GenomeManager

        digest = _build_on_disk_store(tmp_path / "vgp", fixtures_path / "rCRSd.fa", "rCRSd")
        rg = _server_rgc()
        rg.store.add("vgp", str(tmp_path / "vgp"), StoreType.on_disk, priority=10)
        srow = rg.store.get("vgp")
        priorities = {"vgp": 10}

        _sync_one_store(rg, srow, priorities, 100)
        genome = rg.genome.get(digest)
        assert genome.species_name is None
        assert genome.taxon_id is None

        fhr = {
            "genome": "Homo sapiens",
            "commonName": "human",
            "taxon": {"name": "Homo sapiens", "uri": "https://identifiers.org/taxonomy:9606"},
            "assemblySource": "NCBI",
            "accessionID": {"name": "GCA_000001405.15"},
            "assemblyLevel": "chromosome",
        }
        store = RefgetStore.on_disk(str(tmp_path / "vgp"))
        store.set_fhr_metadata(digest, GenomeManager._fhr_metadata_from_dict(fhr))

        _sync_one_store(rg, srow, priorities, 100)
        genome = rg.genome.get(digest)
        assert genome.species_name == "Homo sapiens"
        assert genome.common_name == "human"
        assert genome.taxon_id == 9606
        assert genome.assembly_source == "NCBI"
        assert genome.assembly_accession == "GCA_000001405.15"
        assert genome.assembly_level == "chromosome"


# ---------------------------------------------------------------------------
# Local mode reads both halves of the alias space
# ---------------------------------------------------------------------------

#: A digest for a genome this node did not build. Only its SQL rows exist here.
FEDERATED_DIGEST = "F" * 32


def _local_rg(tmp_path):
    """A build-node-like instance: local mode, one real on-disk store."""
    rg = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    rg.init(genome_folder=tmp_path / "genomes", genome_stage_folder=tmp_path / "stage")
    return rg


def _register_federated(rg, alias, digest=FEDERATED_DIGEST, store="vgp"):
    """Register a genome and alias the way ``store sync`` does: SQL only."""
    from refgenie.exceptions import MissingGenomeError

    if not any(s.name == store for s in rg.store.list_all()):
        rg.store.add(store, f"url://{store}", StoreType.remote, priority=20)
    try:
        rg.genome.get(digest)
    except MissingGenomeError:
        rg.genome.add(digest, "federated genome", [], store_name=store)
    priorities = {s.name: s.priority for s in rg.store.list_all()}
    AliasManager(rg.database_engine).federated_sync(alias, digest, store, priorities)
    return digest


class TestFederatedAliasReads:
    def test_resolve_finds_federated_only_name(self, tmp_path):
        rg = _local_rg(tmp_path)
        _register_federated(rg, "rAllMis1")
        assert rg.alias.resolve("rAllMis1") == FEDERATED_DIGEST

    def test_local_alias_wins_over_federated(self, tmp_path, fixtures_path):
        rg = _local_rg(tmp_path)
        local_digest, _ = rg.genome.initialize_genome(
            fasta_file_path=fixtures_path / "rCRSd.fa",
            alias_names=["rCRSd"],
            description="built here",
        )
        _register_federated(rg, "rCRSd")

        assert local_digest != FEDERATED_DIGEST
        # The name belongs to the genome this node built.
        assert rg.alias.resolve("rCRSd") == local_digest
        # The federated genome is still reachable by digest and store::name.
        assert rg.genome.get(FEDERATED_DIGEST).digest == FEDERATED_DIGEST
        assert rg.alias.resolve("vgp::rCRSd") == FEDERATED_DIGEST

    def test_list_all_unions_without_duplicates(self, tmp_path, fixtures_path):
        rg = _local_rg(tmp_path)
        local_digest, _ = rg.genome.initialize_genome(
            fasta_file_path=fixtures_path / "rCRSd.fa",
            alias_names=["rCRSd"],
            description="built here",
        )
        _register_federated(rg, "rCRSd")
        _register_federated(rg, "rAllMis1")

        pairs = [(a.name, a.genome_digest) for a in rg.alias.list_all()]
        # One row per name: the contested 'rCRSd' appears once, pointing at the
        # genome this node built.
        assert len(pairs) == len({name for name, _ in pairs})
        assert set(pairs) == {
            ("rCRSd", local_digest),
            ("rAllMis1", FEDERATED_DIGEST),
        }

    def test_get_for_genome_unions(self, tmp_path):
        rg = _local_rg(tmp_path)
        _register_federated(rg, "rAllMis1")
        assert rg.alias.get_for_genome(FEDERATED_DIGEST) == ["rAllMis1"]

    def test_exists_sees_federated_name(self, tmp_path):
        rg = _local_rg(tmp_path)
        _register_federated(rg, "rAllMis1")
        assert rg.alias.exists("rAllMis1")
        assert not rg.alias.exists("nobody")

    def test_add_refuses_a_federated_name(self, tmp_path, fixtures_path):
        rg = _local_rg(tmp_path)
        _register_federated(rg, "rAllMis1")
        local_digest, _ = rg.genome.initialize_genome(
            fasta_file_path=fixtures_path / "rCRSd.fa",
            alias_names=["rCRSd"],
            description="built here",
        )
        with pytest.raises(FederatedAliasError):
            rg.alias.add("rAllMis1", local_digest)
        # The federated genome keeps its name.
        assert rg.alias.resolve("rAllMis1") == FEDERATED_DIGEST

    def test_add_is_idempotent_for_the_same_genome(self, tmp_path):
        rg = _local_rg(tmp_path)
        _register_federated(rg, "rAllMis1")
        rg.alias.add("rAllMis1", FEDERATED_DIGEST)
        assert rg.alias.resolve("rAllMis1") == FEDERATED_DIGEST

    def test_no_stores_behaves_like_the_store_alone(self, tmp_path, fixtures_path):
        """With an empty registry the SQL half contributes nothing."""
        rg = _local_rg(tmp_path)
        local_digest, _ = rg.genome.initialize_genome(
            fasta_file_path=fixtures_path / "rCRSd.fa",
            alias_names=["rCRSd"],
            description="built here",
        )
        assert [(a.name, a.genome_digest) for a in rg.alias.list_all()] == [
            ("rCRSd", local_digest)
        ]
        assert rg.alias.resolve("rCRSd") == local_digest
        with pytest.raises(MissingAliasError):
            rg.alias.resolve("nobody")

    def test_reload_store_router_invalidates_the_alias_cache(self, tmp_path, fixtures_path):
        rg = _local_rg(tmp_path)
        rg.genome.initialize_genome(
            fasta_file_path=fixtures_path / "rCRSd.fa",
            alias_names=["rCRSd"],
            description="built here",
        )
        assert len(rg.alias.list_all()) == 1
        # Remove the alias behind the manager's back, as another process would.
        rg.refget_store.remove_collection_alias("refgenie", "rCRSd")
        rg.reload_store_router()
        assert rg.alias.list_all() == []
