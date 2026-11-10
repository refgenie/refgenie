"""Federated multi-store serving: store registry, alias policy, sync, and serving
a locally built database.

These build real on-disk RefgetStores and a real catalog, so they are
``component`` tier.
"""

import json
import logging
import socket

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
from refgenie.managers.store import LOCAL_STORE_NAME
from refgenie.models import GenomeAlias

from tests.support.factories import build_on_disk_store, make_engine, seed_species
from tests.support.server import needs_server

pytestmark = pytest.mark.component


def _server_rgc():
    rg = Refgenie(database_engine=make_engine(), suppress_migrations=True, server_mode=True)
    rg.database.create_tables()
    return rg


def _write_species(store_dir, species):
    """Write a species.json at an on-disk store root."""
    (store_dir / "species.json").write_text(json.dumps({"species": species}))


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

    def test_add_reserved_local_name_raises(self):
        rg = _server_rgc()
        with pytest.raises(ValueError, match="reserved"):
            rg.store.add(LOCAL_STORE_NAME, "url://l", StoreType.remote)
        assert list(rg.store.list_all()) == []

    def test_env_declaration_skips_reserved_local_name(self, caplog):
        rg = _server_rgc()
        rg.store.apply_declaration(
            [{"name": LOCAL_STORE_NAME, "url": "url://l"}, {"name": "vgp", "url": "url://v"}],
        )
        assert [s.name for s in rg.store.list_all()] == ["vgp"]
        assert "reserved" in caplog.text

    def test_apply_declaration_reports_changes(self):
        rg = _server_rgc()
        rg.store.add("jungle", "url://j", StoreType.remote, priority=10)
        rg.store.add("vgp", "url://v", StoreType.remote, priority=20)
        changes = rg.store.apply_declaration(
            [
                {"name": "jungle", "url": "url://j", "priority": 10},
                {"name": "vgp", "url": "url://v2", "priority": 20},
                {"name": "igenomes", "url": "url://i", "priority": 40},
            ]
        )
        assert changes == {"added": ["igenomes"], "updated": ["vgp"], "disabled": []}
        assert rg.store.get("vgp").url == "url://v2"

    def test_apply_declaration_keeps_undeclared_stores_by_default(self):
        rg = _server_rgc()
        rg.store.add("jungle", "url://j", StoreType.remote, priority=10)
        rg.store.apply_declaration([{"name": "vgp", "url": "url://v"}])
        assert rg.store.get("jungle").enabled

    def test_apply_declaration_disables_missing(self):
        """Undeclared stores are disabled, never deleted: genomes keep their FK."""
        rg = _server_rgc()
        rg.store.add("jungle", "url://j", StoreType.remote, priority=10)
        rg.store.add("old", "url://o", StoreType.remote, priority=50)
        rg.genome.add("D" * 32, "owned by old", [], store_name="old")

        changes = rg.store.apply_declaration(
            [{"name": "jungle", "url": "url://j", "priority": 10}], disable_missing=True
        )

        assert changes == {"added": [], "updated": [], "disabled": ["old"]}
        assert [s.name for s in rg.store.list_all()] == ["jungle", "old"]
        assert [s.name for s in rg.store.enabled_stores()] == ["jungle"]
        assert rg.genome.get("D" * 32).store_name == "old"
        # Applying the same file again changes nothing.
        again = rg.store.apply_declaration(
            [{"name": "jungle", "url": "url://j", "priority": 10}], disable_missing=True
        )
        assert again == {"added": [], "updated": [], "disabled": []}

    def test_apply_declaration_reenables_a_declared_store(self):
        rg = _server_rgc()
        rg.store.add("vgp", "url://v", StoreType.remote, priority=20, enabled=False)
        changes = rg.store.apply_declaration(
            [{"name": "vgp", "url": "url://v", "priority": 20}], disable_missing=True
        )
        assert changes["updated"] == ["vgp"]
        assert rg.store.get("vgp").enabled

    def test_store_apply_cli_reads_the_file(self, tmp_path, capsys):
        from refgenie.cli.commands.store import StoreApplyModel, handle_store_apply

        rg = _server_rgc()
        rg.store.add("old", "url://o", StoreType.remote, priority=50)
        declaration = tmp_path / "stores.json"
        declaration.write_text(
            json.dumps([{"name": "jungle", "url": "url://j", "priority": 10, "type": "remote"}])
        )

        handle_store_apply(StoreApplyModel(file=declaration), rg)

        assert [s.name for s in rg.store.enabled_stores()] == ["jungle"]
        out = capsys.readouterr().out
        assert "added" in out and "disabled" in out

    def test_store_apply_cli_rejects_a_non_array(self, tmp_path):
        from refgenie.cli.commands.store import StoreApplyModel, handle_store_apply

        declaration = tmp_path / "stores.json"
        declaration.write_text(json.dumps({"name": "jungle", "url": "url://j"}))
        with pytest.raises(SystemExit):
            handle_store_apply(StoreApplyModel(file=declaration), _server_rgc())

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
        jungle_shared = build_on_disk_store(
            tmp_path / "jungle", fixtures_path / "rCRSd.fa", "shared"
        )
        vgp_shared = build_on_disk_store(tmp_path / "vgp", fixtures_path / "t7.fa", "shared")

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

    def test_sync_makes_no_network_call(self, tmp_path, fixtures_path, monkeypatch):
        """Syncing an on-disk store opens no connection: nothing, NCBI included,
        is fetched. Its species.json is read from disk."""
        from refgenie.cli.commands.store import StoreSyncModel, handle_store_sync

        digest = build_on_disk_store(tmp_path / "jungle", fixtures_path / "rCRSd.fa", "rCRSd")
        _write_species(tmp_path / "jungle", [{"taxon_id": 10116, "scientific_name": "Rattus"}])
        rg = _server_rgc()
        rg.store.add("jungle", str(tmp_path / "jungle"), StoreType.on_disk, priority=10)

        def no_network(*args, **kwargs):
            raise AssertionError(f"store sync opened a connection: {args}")

        monkeypatch.setattr(socket, "create_connection", no_network)
        monkeypatch.setattr(socket.socket, "connect", no_network)
        handle_store_sync(StoreSyncModel(), rg)

        assert rg.genome.get(digest).store_name == "jungle"
        assert rg.species.get(10116).scientific_name == "Rattus"

    def test_single_store_sync_reloads_species_from_every_enabled_store(
        self, tmp_path, fixtures_path
    ):
        from refgenie.cli.commands.store import StoreSyncModel, handle_store_sync

        build_on_disk_store(tmp_path / "jungle", fixtures_path / "rCRSd.fa", "rCRSd")
        build_on_disk_store(tmp_path / "vgp", fixtures_path / "t7.fa", "t7")
        _write_species(
            tmp_path / "jungle",
            [{"taxon_id": 1, "scientific_name": "One", "common_name": "from jungle"}],
        )
        _write_species(
            tmp_path / "vgp",
            [
                {"taxon_id": 1, "scientific_name": "One", "common_name": "from vgp"},
                {"taxon_id": 2, "scientific_name": "Two"},
            ],
        )
        rg = _server_rgc()
        rg.store.add("jungle", str(tmp_path / "jungle"), StoreType.on_disk, priority=10)
        rg.store.add("vgp", str(tmp_path / "vgp"), StoreType.on_disk, priority=20)

        handle_store_sync(StoreSyncModel(name="jungle"), rg)

        assert rg.species.get(2) is not None
        assert rg.species.get(1).common_name == "from jungle"

    def test_identical_genome_across_stores_dedups(self, tmp_path, fixtures_path):
        """The same FASTA in two stores yields one Genome row (digest dedup)."""
        from refgenie.cli.commands.store import _sync_one_store

        digest_a = build_on_disk_store(tmp_path / "a", fixtures_path / "rCRSd.fa", "rCRSd")
        digest_b = build_on_disk_store(tmp_path / "b", fixtures_path / "rCRSd.fa", "rCRSd")
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

    def test_disabled_owner_loses_genome(self, tmp_path, fixtures_path):
        """A genome owned by a disabled store moves to an enabled store that offers it."""
        from refgenie.cli.commands.store import StoreSyncModel, handle_store_sync

        digest = build_on_disk_store(tmp_path / "new", fixtures_path / "rCRSd.fa", "rCRSd")
        rg = _server_rgc()
        rg.store.add("old", "url://old", StoreType.remote, priority=10, enabled=False)
        rg.store.add("new", str(tmp_path / "new"), StoreType.on_disk, priority=20)
        rg.genome.add(digest=digest, description="", alias_names=[], store_name="old")

        handle_store_sync(StoreSyncModel(), rg)

        assert rg.genome.get(digest).store_name == "new"

    def test_resync_backfills_fhr_columns(self, tmp_path, fixtures_path):
        """A sidecar added to the store after the first sync lands on re-sync."""
        from refget.store import RefgetStore

        from refgenie.cli.commands.store import _sync_one_store
        from refgenie.managers.genome import GenomeManager

        digest = build_on_disk_store(tmp_path / "vgp", fixtures_path / "rCRSd.fa", "rCRSd")
        rg = _server_rgc()
        rg.store.add("vgp", str(tmp_path / "vgp"), StoreType.on_disk, priority=10)
        srow = rg.store.get("vgp")
        priorities = {"vgp": 10}

        _sync_one_store(rg, srow, priorities, 100)
        genome = rg.genome.get(digest)
        assert genome.species_name is None
        assert genome.taxon_id is None

        fhr = {
            "genome": "rCRSd",
            "genomeSynonym": ["rCRS"],
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
        assert genome.synonyms == ["rCRS"]
        assert genome.taxon_id == 9606
        assert genome.assembly_source == "NCBI"
        assert genome.assembly_accession == "GCA_000001405.15"
        assert genome.assembly_level == "chromosome"
        seed_species(rg, 9606, "Homo sapiens", "human")
        assert rg.genome.detail(digest).common_name == "human"

        del fhr["genomeSynonym"]
        store.set_fhr_metadata(digest, GenomeManager._fhr_metadata_from_dict(fhr))
        _sync_one_store(rg, srow, priorities, 100)
        assert rg.genome.get(digest).synonyms is None


# ---------------------------------------------------------------------------
# Comparing genomes that live in different stores
# ---------------------------------------------------------------------------


def _write_fasta(path, records):
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records))
    return path


def _synced_rgc(stores):
    """A server-mode refgenie with the given (name, path) on-disk stores synced."""
    from refgenie.cli.commands.store import _sync_one_store

    rg = _server_rgc()
    for priority, (name, path) in enumerate(stores):
        rg.store.add(name, str(path), StoreType.on_disk, priority=priority)
    priorities = {s.name: s.priority for s in rg.store.list_all()}
    for srow in rg.store.enabled_stores():
        _sync_one_store(rg, srow, priorities, 100)
    return rg


class TestCompareAcrossStores:
    @pytest.fixture
    def fastas(self, tmp_path):
        # Two shared sequences in the same order, one unique to each genome.
        a = _write_fasta(
            tmp_path / "a.fa", [("chr1", "ACGTACGT"), ("chr2", "GGGGCCCC"), ("chr3", "TTTTAA")]
        )
        b = _write_fasta(
            tmp_path / "b.fa", [("chr1", "ACGTACGT"), ("chr2", "GGGGCCCC"), ("chr4", "CA")]
        )
        return a, b

    def test_cross_store_compare_returns_the_seqcol_comparison(self, tmp_path, fastas):
        digest_a = build_on_disk_store(tmp_path / "s1", fastas[0], "a")
        digest_b = build_on_disk_store(tmp_path / "s2", fastas[1], "b")
        rg = _synced_rgc([("s1", tmp_path / "s1"), ("s2", tmp_path / "s2")])
        assert rg.genome.get(digest_a).store_name == "s1"
        assert rg.genome.get(digest_b).store_name == "s2"

        result = rg.genome.compare(digest_a, digest_b)

        assert result["digests"] == {"a": digest_a, "b": digest_b}
        assert set(result["attributes"]) == {"a_only", "b_only", "a_and_b"}
        assert {"names", "lengths", "sequences", "name_length_pairs"} <= set(
            result["attributes"]["a_and_b"]
        )
        elements = result["array_elements"]
        assert set(elements) == {"a_count", "b_count", "a_and_b_count", "a_and_b_same_order"}
        assert elements["a_count"]["names"] == 3
        assert elements["a_and_b_count"]["sequences"] == 2
        assert elements["a_and_b_same_order"]["names"] is True

    def test_same_store_and_cross_store_agree(self, tmp_path, fastas):
        build_on_disk_store(tmp_path / "both", fastas[0], "a")
        build_on_disk_store(tmp_path / "both", fastas[1], "b")
        digest_a = build_on_disk_store(tmp_path / "s1", fastas[0], "a")
        digest_b = build_on_disk_store(tmp_path / "s2", fastas[1], "b")

        same = _synced_rgc([("both", tmp_path / "both")])
        cross = _synced_rgc([("s1", tmp_path / "s1"), ("s2", tmp_path / "s2")])

        assert same.genome.compare(digest_a, digest_b) == cross.genome.compare(digest_a, digest_b)


# ---------------------------------------------------------------------------
# Local mode reads both halves of the alias space
# ---------------------------------------------------------------------------

#: A digest for a genome this node did not build. Only its SQL rows exist here.
FEDERATED_DIGEST = "F" * 32


def _local_rg(tmp_path):
    """A build-node-like instance: local mode, one real on-disk store."""
    rg = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    rg.database.init(genome_folder=tmp_path / "genomes", genome_stage_folder=tmp_path / "stage")
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
        assert [(a.name, a.genome_digest) for a in rg.alias.list_all()] == [("rCRSd", local_digest)]
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


# ---------------------------------------------------------------------------
# Serving a locally built database (`refgenie serve` runs in server mode)
# ---------------------------------------------------------------------------


def _serve_mode_twin(rg):
    """The server-mode instance `refgenie serve` would open over ``rg``'s database."""
    return Refgenie(database_engine=rg.database_engine, suppress_migrations=True, server_mode=True)


class TestServingALocalBuild:
    """A database built on this node keeps its aliases in its own RefgetStore.

    Server mode must still see them: that store joins the federation and its
    aliases are read before the SQL table's.
    """

    pytestmark = needs_server

    def test_aliases_come_from_the_nodes_own_store(self, tmp_path, fixtures_path):
        built = _local_rg(tmp_path)
        digest, _ = built.genome.initialize_genome(
            fasta_file_path=fixtures_path / "rCRSd.fa",
            alias_names=["rCRSd"],
            description="built here",
        )
        served = _serve_mode_twin(built)
        assert served.alias.resolve("rCRSd") == digest
        assert served.store_router.names[0] == "local"

    def test_local_name_wins_over_a_federated_one(self, tmp_path, fixtures_path):
        from refgenie.cli.commands.store import _sync_one_store

        built = _local_rg(tmp_path)
        local_digest, _ = built.genome.initialize_genome(
            fasta_file_path=fixtures_path / "rCRSd.fa",
            alias_names=["rCRSd"],
            description="built here",
        )
        # A federated store names a different genome 'rCRSd', and one 'phage'.
        vgp = tmp_path / "vgp"
        t7_digest = build_on_disk_store(vgp, fixtures_path / "t7.fa", "rCRSd")
        build_on_disk_store(vgp, fixtures_path / "t7.fa", "phage")
        built.store.add("vgp", str(vgp), StoreType.on_disk, priority=10)
        _sync_one_store(built, built.store.get("vgp"), {"vgp": 10}, 100)

        served = _serve_mode_twin(built)
        assert served.store_router.names == ["local", "vgp"]
        assert served.alias.resolve("rCRSd") == local_digest
        assert served.alias.resolve("phage") == t7_digest
        assert served.alias.resolve("vgp::rCRSd") == t7_digest

    def test_no_own_store_reads_sql_only(self, tmp_path):
        rg = _local_rg(tmp_path)
        _register_federated(rg, "rAllMis1")
        served = _serve_mode_twin(rg)
        assert not (served.genome_folder / ".refget_store").exists()
        assert isinstance(served.alias, AliasManager)
        assert served.alias.resolve("rAllMis1") == FEDERATED_DIGEST

    def test_served_endpoints_list_and_resolve_the_alias(self, tmp_path, fixtures_path):
        from tests.support.factories import make_built_refgenie
        from tests.support.server import make_server_client

        served = _serve_mode_twin(make_built_refgenie(tmp_path, fixtures_path))
        with make_server_client(served) as client:
            genomes = client.get("/v4/genomes").json()["results"]
            assert [g["aliases"] for g in genomes] == [["rCRSd"]]
            response = client.get("/v4/aliases/rCRSd")
        assert response.status_code == 200
        assert response.json()["digest"] == genomes[0]["digest"]

    def test_client_pulls_by_alias(self, tmp_path, fixtures_path):
        from tests.support.factories import MOCK_SERVER_URL
        from tests.support.server import make_server_client_world, serve_refgenie

        built, client = make_server_client_world(make_engine(), tmp_path, fixtures_path)
        with serve_refgenie(client, _serve_mode_twin(built), MOCK_SERVER_URL):
            pulled = client.transfer.pull(
                genome=GenomeAlias("rCRSd"),
                asset_group_name="fasta",
                force=True,
                force_large=True,
                force_server_urls=[MOCK_SERVER_URL],
            )
        assert pulled is not None
        assert client.alias.resolve("rCRSd") == built.alias.resolve("rCRSd")


# ---------------------------------------------------------------------------
# Server mode: one read-only store per dataset, shared with seqcol
# ---------------------------------------------------------------------------


class TestSharedReadonlyStores:
    """The router's frozen stores are the very objects the seqcol mounts serve."""

    pytestmark = needs_server

    FHR = {
        "genome": "rCRSd",
        "taxon": {"name": "Homo sapiens", "uri": "https://identifiers.org/taxonomy:9606"},
    }

    def _served(self, tmp_path, fixtures_path):
        from refget.store import RefgetStore

        from refgenie.cli.commands.store import _sync_one_store
        from refgenie.managers.genome import GenomeManager

        human = build_on_disk_store(tmp_path / "jungle", fixtures_path / "rCRSd.fa", "rCRSd")
        phage = build_on_disk_store(tmp_path / "vgp", fixtures_path / "t7.fa", "t7")
        RefgetStore.on_disk(str(tmp_path / "jungle")).set_fhr_metadata(
            human, GenomeManager._fhr_metadata_from_dict(self.FHR)
        )
        rg = _server_rgc()
        rg.store.add("jungle", str(tmp_path / "jungle"), StoreType.on_disk, priority=10)
        rg.store.add("vgp", str(tmp_path / "vgp"), StoreType.on_disk, priority=20)
        priorities = {"jungle": 10, "vgp": 20}
        for srow in rg.store.enabled_stores():
            _sync_one_store(rg, srow, priorities, 100)
        return rg, human, phage

    @staticmethod
    def _mounts(app):
        return {m.url: m.app for m in app.state.seqcol.mounts}

    def test_seqcol_mounts_share_the_routers_readonly_stores(
        self, tmp_path, fixtures_path, monkeypatch
    ):
        from refget.middleware import StoreFreshnessMiddleware
        from refget.store import ReadonlyRefgetStore

        from tests.support.server import make_server_client

        rg, human, phage = self._served(tmp_path, fixtures_path)
        # Raising (the default): the endpoint-hit middleware runs after every
        # /seqcol response, and must not error inside a mounted sub-app.
        with make_server_client(rg) as client:
            mounts = self._mounts(client.app)
            router = rg.store_router
            for path, name in (("/seqcol", "jungle"), ("/seqcol/vgp", "vgp")):
                store = router.get_store(name)
                assert isinstance(store, ReadonlyRefgetStore)
                assert mounts[path].state.backend._store is store
                assert not any(
                    m.cls is StoreFreshnessMiddleware for m in mounts[path].user_middleware
                )
            # Exactly one store object per dataset.
            assert len({id(s) for _, s in router.iter_stores()}) == len(router.names) == 2

            # Nothing opens or loads a store once startup is done.
            def _no_io(*a, **k):
                raise AssertionError("store opened or loaded after startup")

            monkeypatch.setattr("refgenie.core.store_router.make_store", _no_io)
            monkeypatch.setattr("refget.store.RefgetStore.on_disk", _no_io)
            monkeypatch.setattr("refget.store.RefgetStore.open_remote", _no_io)
            monkeypatch.setattr("refget.seqcolapi.app.prepare_store", _no_io)

            detail = client.get(f"/v4/genomes/{human}")
            assert detail.status_code == 200
            assert detail.json()["fhr"]["genome"] == "rCRSd"
            assert detail.json()["species_name"] == "Homo sapiens"
            alias = client.get("/v4/aliases/rCRSd")
            assert alias.status_code == 200
            assert alias.json()["fhr"]["genome"] == "rCRSd"
            assert client.get(f"/v4/genomes/{phage}").json()["fhr"] is None

            coll = client.get(f"/seqcol/collection/{human}?level=2")
            assert coll.status_code == 200
            assert coll.json()["names"]
            coll = client.get(f"/seqcol/vgp/collection/{phage}?level=2")
            assert coll.status_code == 200
            cmp = client.get(f"/seqcol/comparison/{human}/{human}")
            assert cmp.status_code == 200
            assert cmp.json()["digests"]["a"] == human

            # Mounted seqcol hits are counted under the host app's collector.
            collector = client.app.state.endpoint_hits_collector
            assert [s.total_hits for s in collector.get_stats("collection").values()] == [1, 1]
            assert [s.total_hits for s in collector.get_stats("compare_2_digests").values()] == [1]
