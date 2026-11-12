"""The publish catalog (``refgenie.publish_catalog``): export, import and prune."""

import pathlib

import pytest
from sqlalchemy import create_engine
from sqlmodel import SQLModel, select

from refgenie import Refgenie
from refgenie.const import TARGET_ALEMBIC_VERSION
from refgenie.db.tables import (
    AlembicVersion,
    Alias,
    Asset,
    Genome,
    Remote,
    RemoteAssetLink,
    RemoteType,
    StagedAsset,
    Store,
    StoreType,
)
from refgenie.managers.alias import AliasManager, _StoreAliases
from refgenie.publish_catalog import (
    export_publish_catalog,
    import_publish_catalog,
    prune_publish_catalog,
)
from tests.support.factories import make_engine, register_fasta, seed_alias_row
from tests.support.server import make_server_client, requires_server

pytestmark = pytest.mark.component

HTTPS_PREFIX = "https://example.org/assets"

#: A genome this node did not build. Only its SQL rows exist here.
FEDERATED_DIGEST = "F" * 32


def artifact_engine(path):
    return create_engine(f"sqlite:///{path}")


def _build_node(root, fixtures_path):
    """A local-mode catalog holding one locally built and one federated genome."""
    rg = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    rg.database.init(genome_folder=root / "genomes", genome_stage_folder=root / "stage")
    local_digest, _ = rg.genome.initialize_genome(
        fasta_file_path=fixtures_path / "rCRSd.fa",
        alias_names=["rCRSd"],
        description="built here",
    )

    rg.store.add("vgp", "url://vgp", StoreType.remote, priority=20)
    rg.genome.add(FEDERATED_DIGEST, "federated genome", [], store_name="vgp")
    AliasManager(rg.database_engine).federated_sync(
        "rAllMis1", FEDERATED_DIGEST, "vgp", {"vgp": 20}
    )
    return rg, local_digest


@pytest.fixture
def build_node(tmp_path, fixtures_path):
    return _build_node(tmp_path, fixtures_path)


def _new_server(root):
    server = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    server.database.init(genome_folder=root / "server_genomes")
    return server


def _exported_aliases(path):
    from sqlalchemy import create_engine

    engine = create_engine(f"sqlite:///{path}")
    with engine.connect() as conn:
        return {(a.name, a.genome_digest) for a in conn.execute(select(Alias.__table__))}


class TestAliasCoverage:
    def test_export_carries_both_alias_halves(self, build_node, tmp_path):
        rg, local_digest = build_node
        dest = tmp_path / "publish.sqlite"

        summary = export_publish_catalog(rg, dest, HTTPS_PREFIX)

        assert summary["genome"] == 2
        assert summary["alias"] == 2
        assert _exported_aliases(dest) == {
            ("rCRSd", local_digest),
            ("rAllMis1", FEDERATED_DIGEST),
        }

    def test_export_refuses_to_drop_federated_names(self, build_node, tmp_path):
        """The export fails when the alias manager reads only the store half,
        which no table row count would reveal."""
        rg, _local_digest = build_node
        store_only = _StoreAliases(lambda: rg.refget_store, lambda: rg.genome_folder)
        rg._alias_manager = store_only

        with pytest.raises(RuntimeError, match="alias rows the catalog holds"):
            export_publish_catalog(rg, tmp_path / "publish.sqlite", HTTPS_PREFIX)


def _server_aliases(server):
    with server.database_engine.connect() as conn:
        return {(a.name, a.genome_digest) for a in conn.execute(select(Alias.__table__))}


class TestAliasImport:
    """An import rewrites the alias set of each genome the catalog carries."""

    def test_import_replaces_carried_aliases_and_skips_slash_names(self, build_node, tmp_path):
        rg, local_digest = build_node
        dest = tmp_path / "publish.sqlite"
        export_publish_catalog(rg, dest, HTTPS_PREFIX, **CATALOG_STORE)
        server = _new_server(tmp_path)
        server.genome.add("S" * 32, "server only", [], store_name=None)
        seed_alias_row(server, "kept", "S" * 32)
        seed_alias_row(server, "rAllMis1", "S" * 32)  # the catalog assigns it to F
        import_publish_catalog(server.database_engine, str(dest))
        seed_alias_row(server, "renamed-upstream", FEDERATED_DIGEST)
        with artifact_engine(dest).begin() as conn:
            conn.execute(
                Alias.__table__.insert().values(
                    name="Homo_sapiens/NCBI/x", genome_digest=FEDERATED_DIGEST
                )
            )

        summary = import_publish_catalog(server.database_engine, str(dest))

        assert (summary["alias"], summary["alias_skipped"]) == (2, 1)
        assert _server_aliases(server) == {
            ("rCRSd", local_digest),
            ("rAllMis1", FEDERATED_DIGEST),
            ("kept", "S" * 32),
        }
        assert _genomes(server) == {local_digest, FEDERATED_DIGEST, "S" * 32}

    def test_import_refuses_an_artifact_short_of_aliases(self, build_node, tmp_path):
        rg, _local_digest = build_node
        dest = tmp_path / "publish.sqlite"
        export_publish_catalog(rg, dest, HTTPS_PREFIX, **CATALOG_STORE)
        server = _new_server(tmp_path)
        import_publish_catalog(server.database_engine, str(dest))
        for name in ("old-1", "old-2", "old-3"):
            seed_alias_row(server, name, FEDERATED_DIGEST)
        before = _server_aliases(server)

        with pytest.raises(ValueError, match="2 aliases for its genomes, down from 5"):
            import_publish_catalog(server.database_engine, str(dest))
        assert _server_aliases(server) == before


class TestPluginSettings:
    """Plugin settings hold build-node paths; they never reach a published catalog."""

    def test_export_strips_and_import_preserves_plugin_settings(self, build_node, tmp_path):
        from refgenie.publish_catalog import import_publish_catalog

        rg, _local_digest = build_node
        rg.plugins.set("nfcore", config_path="/build/node/nf.config")
        dest = tmp_path / "publish.sqlite"
        export_publish_catalog(rg, dest, HTTPS_PREFIX, **CATALOG_STORE)

        exported = Refgenie(
            database_engine=artifact_engine(dest), suppress_migrations=True
        ).plugins.all_settings()
        assert exported == {}

        server = _new_server(tmp_path)
        server.plugins.set("own", key="kept")
        import_publish_catalog(server.database_engine, str(dest))
        assert server.plugins.all_settings() == {"own": {"key": "kept"}}


def _exported_stores(path):
    with artifact_engine(path).connect() as conn:
        return {s.name: s for s in conn.execute(select(Store.__table__))}


def _exported_genome_stores(path):
    with artifact_engine(path).connect() as conn:
        return {g.digest: g.store_name for g in conn.execute(select(Genome.__table__))}


CATALOG_STORE = {
    "catalog_store_name": "refgenie-main",
    "catalog_store_url": "https://example.org/refget-store/refgenie-main/",
}


class TestStores:
    """The store list travels in the catalog, so a server learns its stores from it."""

    def test_export_publishes_remote_stores_only(self, build_node, tmp_path):
        rg, local_digest = build_node
        rg.store.add("scratch", str(tmp_path / "scratch"), StoreType.on_disk, priority=5)
        dest = tmp_path / "publish.sqlite"

        summary = export_publish_catalog(rg, dest, HTTPS_PREFIX)

        assert summary["store"] == 1
        assert set(_exported_stores(dest)) == {"vgp"}
        assert _exported_genome_stores(dest) == {local_digest: None, FEDERATED_DIGEST: "vgp"}

    def test_catalog_store_owns_the_nodes_own_genomes(self, build_node, tmp_path):
        rg, local_digest = build_node
        dest = tmp_path / "publish.sqlite"

        export_publish_catalog(rg, dest, HTTPS_PREFIX, catalog_store_priority=15, **CATALOG_STORE)

        stores = _exported_stores(dest)
        assert set(stores) == {"vgp", "refgenie-main"}
        main = stores["refgenie-main"]
        assert (main.url, main.type, main.priority, main.enabled) == (
            CATALOG_STORE["catalog_store_url"],
            StoreType.remote,
            15,
            True,
        )
        assert _exported_genome_stores(dest) == {
            local_digest: "refgenie-main",
            FEDERATED_DIGEST: "vgp",
        }
        # The build catalog itself is not rewritten.
        assert rg.genome.get(local_digest).store_name is None

    def test_export_refuses_a_genome_in_an_unpublished_store(self, build_node, tmp_path):
        """A genome in an on-disk store would break the server's store foreign key."""
        rg, _ = build_node
        rg.store.add("scratch", str(tmp_path / "scratch"), StoreType.on_disk, priority=5)
        rg.genome.add("S" * 32, "in a build-node path", [], store_name="scratch")
        dest = tmp_path / "publish.sqlite"
        dest.write_text("previous artifact")

        with pytest.raises(RuntimeError, match="'scratch' \\(1 genome\\)"):
            export_publish_catalog(rg, dest, HTTPS_PREFIX, **CATALOG_STORE)
        assert dest.read_text() == "previous artifact"

    def test_export_needs_both_catalog_store_name_and_url(self, build_node, tmp_path):
        rg, _ = build_node
        with pytest.raises(ValueError, match="both a name and a url"):
            export_publish_catalog(
                rg, tmp_path / "p.sqlite", HTTPS_PREFIX, catalog_store_name="refgenie-main"
            )

    def test_import_upserts_stores_and_disables_absent_ones(self, build_node, tmp_path):
        rg, local_digest = build_node
        dest = tmp_path / "publish.sqlite"
        export_publish_catalog(rg, dest, HTTPS_PREFIX, **CATALOG_STORE)
        server = _new_server(tmp_path)
        server.store.add("vgp", "url://stale", StoreType.remote, priority=99)
        server.store.add("retired", "url://r", StoreType.remote, priority=60)

        summary = import_publish_catalog(server.database_engine, str(dest))

        assert summary["store"] == 2
        assert summary["store_disabled"] == 1
        assert [s.name for s in server.store.enabled_stores()] == ["vgp", "refgenie-main"]
        assert server.store.get("vgp").url == "url://vgp"
        assert server.store.get("vgp").priority == 20
        assert not server.store.get("retired").enabled
        assert server.genome.get(local_digest).store_name == "refgenie-main"

    def test_import_of_an_artifact_without_stores_leaves_stores_alone(self, build_node, tmp_path):
        rg, _ = build_node
        dest = tmp_path / "publish.sqlite"
        export_publish_catalog(rg, dest, HTTPS_PREFIX)
        with artifact_engine(dest).begin() as conn:
            conn.execute(Genome.__table__.update().values(store_name="jungle"))
            conn.execute(Store.__table__.delete())
        server = _new_server(tmp_path)
        server.store.add("jungle", "url://j", StoreType.remote, priority=10)

        import_publish_catalog(server.database_engine, str(dest))

        assert [s.name for s in server.store.enabled_stores()] == ["jungle"]

    @pytest.mark.parametrize("store_name", [None, "plantref"])
    def test_import_refuses_genomes_without_a_listed_store(self, build_node, tmp_path, store_name):
        from refgenie.db.publish import CatalogStoreError

        rg, _ = build_node
        dest = tmp_path / "publish.sqlite"
        export_publish_catalog(rg, dest, HTTPS_PREFIX, **CATALOG_STORE)
        with artifact_engine(dest).begin() as conn:
            conn.execute(Genome.__table__.update().values(store_name=store_name))
        server = _new_server(tmp_path)

        with pytest.raises(CatalogStoreError, match="could not serve"):
            import_publish_catalog(server.database_engine, str(dest))
        assert list(server.genome.list_all()) == []
        assert list(server.store.enabled_stores()) == []


class TestSpecies:
    """Species tables travel in the catalog; an import rewrites each taxon it carries."""

    @staticmethod
    def _entry(taxon_id, name, synonyms=None):
        from refgenie.utils.species_file import SpeciesEntry

        return SpeciesEntry(taxon_id=taxon_id, scientific_name=name, synonyms=synonyms)

    def _publish(self, rg, server, tmp_path, name):
        from refgenie.publish_catalog import import_publish_catalog

        dest = tmp_path / name
        summary = export_publish_catalog(rg, dest, HTTPS_PREFIX, **CATALOG_STORE)
        import_publish_catalog(server.database_engine, str(dest))
        return summary

    def test_publish_carries_species(self, build_node, tmp_path):
        rg, _ = build_node
        rg.species.replace([self._entry(10116, "Rattus norvegicus", ["brown rat", "rat"])])
        server = _new_server(tmp_path)

        summary = self._publish(rg, server, tmp_path, "publish.sqlite")

        assert summary["species"] == 1
        assert summary["speciesname"] == 3
        assert server.species.get(10116).synonyms == ["brown rat", "rat"]

    def test_reimport_rewrites_species_per_taxon(self, build_node, tmp_path):
        rg, _ = build_node
        rat = self._entry(10116, "Rattus norvegicus", ["brown rat", "rat"])
        rg.species.replace([rat, self._entry(4577, "Zea mays", ["corn"])])
        server = _new_server(tmp_path)
        self._publish(rg, server, tmp_path, "first.sqlite")
        server.species.replace([self._entry(9999, "Server only")])

        rg.species.replace([rat.model_copy(update={"synonyms": ["brown rat"]})])
        self._publish(rg, server, tmp_path, "second.sqlite")

        assert server.species.get(10116).synonyms == ["brown rat"]
        assert server.species.get(4577).synonyms == ["corn"]
        assert server.species.get(9999) is not None


def _add_group(path, digest, group_id, asset_digest):
    from datetime import datetime, timezone

    from refgenie.db.tables import Asset, AssetClass, AssetGroup, AssetName

    now = datetime.now(timezone.utc)
    engine = artifact_engine(path)
    with engine.begin() as conn:
        if conn.execute(select(AssetClass.__table__)).first() is None:
            conn.execute(
                AssetClass.__table__.insert().values(
                    id=1, name="fasta", version="0.0.1", serving_modes=["archive"]
                )
            )
        conn.execute(AssetGroup.__table__.delete())
        conn.execute(
            AssetGroup.__table__.insert().values(
                id=group_id,
                name="bowtie2_index",
                asset_class_id=1,
                genome_digest=digest,
                created_at=now,
                updated_at=now,
            )
        )
        conn.execute(
            Asset.__table__.insert().values(
                digest=asset_digest, name="bowtie2_index", asset_group_id=group_id
            )
        )
        conn.execute(
            AssetName.__table__.insert().values(
                name="default",
                asset_group_id=group_id,
                asset_digest=asset_digest,
                is_default=True,
            )
        )


@pytest.fixture(scope="module")
def artifacts(tmp_path_factory, fixtures_path):
    """Exports of one node, each with a hand-added bowtie2_index group.

    ``first`` holds group 129 with asset a*64. ``second`` is the node after a
    rebuild: group 146 with asset b*64. ``bare`` has no such group, and
    ``empty`` is a node with no genomes. Tests only read these files.
    """
    root = tmp_path_factory.mktemp("artifacts")
    rg, digest = _build_node(root, fixtures_path)
    paths = {name: root / f"{name}.sqlite" for name in ("first", "second", "bare")}
    for path in paths.values():
        export_publish_catalog(rg, path, HTTPS_PREFIX, **CATALOG_STORE)
    _add_group(paths["first"], digest, 129, "a" * 64)
    _add_group(paths["second"], digest, 146, "b" * 64)
    empty = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    empty.database.init(genome_folder=root / "g2", genome_stage_folder=root / "s2")
    paths["empty"] = root / "empty.sqlite"
    export_publish_catalog(empty, paths["empty"], HTTPS_PREFIX, **CATALOG_STORE)
    return digest, paths


def _imported(tmp_path, artifacts, *names):
    server = _new_server(tmp_path)
    for name in names:
        import_publish_catalog(server.database_engine, str(artifacts[1][name]))
    return server


def _assets(server):
    with server.database_engine.connect() as conn:
        return {a.digest for a in conn.execute(select(Asset.__table__))}


def test_reimport_with_new_group_id_keeps_server_id(tmp_path, artifacts):
    """A rebuilt asset's group gets a new id on the build node. The server keeps
    its own id and points the new asset and names at it, rather than fail the
    import on unique_genome_group_name."""
    from refgenie.db.tables import AssetGroup, AssetName

    server = _imported(tmp_path, artifacts, "first", "second")
    with server.database_engine.connect() as conn:
        groups = conn.execute(select(AssetGroup.__table__)).all()
        assets = {a.digest: a.asset_group_id for a in conn.execute(select(Asset.__table__))}
        names = {
            (n.asset_digest, n.asset_group_id) for n in conn.execute(select(AssetName.__table__))
        }
    assert [(g.id, g.name) for g in groups] == [(129, "bowtie2_index")]
    assert assets == {"a" * 64: 129, "b" * 64: 129}
    assert names == {("b" * 64, 129)}


#: A prune result's genome keys when no genome left the artifact.
NOTHING_ABSENT = {"genome": [], "alias": [], "store": [], "kept_genome": []}


def _genomes(server):
    with server.database_engine.connect() as conn:
        return {g.digest for g in conn.execute(select(Genome.__table__))}


def _alias_names(server):
    with server.database_engine.connect() as conn:
        return {a.name for a in conn.execute(select(Alias.__table__))}


def _store_names(server):
    with server.database_engine.connect() as conn:
        return {s.name for s in conn.execute(select(Store.__table__))}


class TestCatalogPrune:
    """Import never deletes genomes or assets; ``catalog-prune`` removes what the
    artifact no longer has."""

    def test_prune_removes_retracted_asset(self, tmp_path, artifacts):
        from refgenie.db.tables import AssetGroup, AssetName

        server = _imported(tmp_path, artifacts, "first", "second")
        assert {"a" * 64, "b" * 64} <= _assets(server)

        removed = prune_publish_catalog(server.database_engine, str(artifacts[1]["second"]))

        assert removed == {"asset": ["a" * 64], "assetgroup": [], **NOTHING_ABSENT}
        assert "a" * 64 not in _assets(server)
        assert "b" * 64 in _assets(server)
        with server.database_engine.connect() as conn:
            names = {n.asset_digest for n in conn.execute(select(AssetName.__table__))}
            groups = [g.name for g in conn.execute(select(AssetGroup.__table__))]
        assert "a" * 64 not in names
        assert groups == ["bowtie2_index"]

    def test_dry_run_deletes_nothing(self, tmp_path, artifacts):
        server = _imported(tmp_path, artifacts, "first", "second")
        before = _assets(server)
        second = str(artifacts[1]["second"])
        removed = prune_publish_catalog(server.database_engine, second, dry_run=True)
        assert removed["asset"] == ["a" * 64]
        assert _assets(server) == before

    def test_prune_drops_group_the_artifact_no_longer_has(self, tmp_path, artifacts):
        from refgenie.db.tables import AssetGroup

        digest, paths = artifacts
        server = _imported(tmp_path, artifacts, "first")

        removed = prune_publish_catalog(server.database_engine, str(paths["bare"]))

        assert removed == {
            "asset": ["a" * 64],
            "assetgroup": [f"{digest}/bowtie2_index"],
            **NOTHING_ABSENT,
        }
        with server.database_engine.connect() as conn:
            assert conn.execute(select(AssetGroup.__table__)).first() is None

    def test_an_artifact_without_genomes_deletes_nothing(self, tmp_path, artifacts):
        server = _imported(tmp_path, artifacts, "first", "second")
        before = _assets(server), _genomes(server)
        removed = prune_publish_catalog(server.database_engine, str(artifacts[1]["empty"]))
        assert removed == {"asset": [], "assetgroup": [], **NOTHING_ABSENT}
        assert (_assets(server), _genomes(server)) == before

    def test_prune_removes_genomes_the_artifact_dropped(self, tmp_path, artifacts):
        """The retired plantref case: a genome with no assets goes, with its alias and
        its now empty disabled store. One that still has assets is kept and listed."""
        import shutil

        digest, paths = artifacts
        server = _imported(tmp_path, artifacts, "first")
        server.store.add("plantref", "url://p", StoreType.remote, priority=50, enabled=False)
        server.genome.add("P" * 32, "retired", [], store_name="plantref")
        seed_alias_row(server, "athal", "P" * 32)
        dropped = shutil.copy(paths["bare"], tmp_path / "dropped.sqlite")
        with artifact_engine(dropped).begin() as conn:
            conn.execute(Alias.__table__.delete().where(Alias.__table__.c.genome_digest == digest))
            conn.execute(Genome.__table__.delete().where(Genome.__table__.c.digest == digest))
        before = _genomes(server), _alias_names(server), _store_names(server)
        expected = {
            "asset": [],
            "assetgroup": [],
            "genome": ["P" * 32],
            "alias": ["athal"],
            "store": ["plantref"],
            "kept_genome": [digest],
        }

        dry = prune_publish_catalog(server.database_engine, str(dropped), dry_run=True)
        assert dry == expected
        assert (_genomes(server), _alias_names(server), _store_names(server)) == before

        assert prune_publish_catalog(server.database_engine, str(dropped)) == expected
        assert _genomes(server) == before[0] - {"P" * 32}
        assert _alias_names(server) == before[1] - {"athal"}
        assert _store_names(server) == before[2] - {"plantref"}


# ---------------------------------------------------------------------------
# A built, staged and pushed asset, exported and imported
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def pushed(tmp_path_factory, fixtures_path):
    """A build node's export: one built, staged asset with one pushed link.

    Yields ``(artifact path, genome digest, asset)``; tests only read the file.
    """
    tmp_path = tmp_path_factory.mktemp("pushed")
    r = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    r.database.init(
        genome_folder=tmp_path / "genomes",
        genome_stage_folder=tmp_path / "stage",
    )
    register_fasta(r, fixtures_path)
    r.genome.initialize_genome(
        fasta_file_path=fixtures_path / "rCRSd.fa",
        alias_names=["rCRSd"],
        description="rCRSd genome",
    )
    r.build.run(
        recipe_name="fasta",
        genome_alias="rCRSd",
        asset_group_name="fasta",
        asset_name="test",
    )
    genome_digest = r.alias.resolve("rCRSd")
    asset = r.asset.get(
        genome_digest=genome_digest,
        asset_group_name="fasta",
        asset_name="test",
    )
    r.stage.create(
        asset=asset,
        genome_folder=r.genome_folder,
        genome_stage_folder=r.genome_stage_folder,
    )
    remote = r.remote.add(
        type=RemoteType.s3,
        prefix="my-bucket/assets",
        name="build-side s3",
        push_command="aws s3 cp {local_path} s3://{prefix}/{relative_path}",
    )
    # archive was pushed; file mode never was
    r.remote.link(ref=remote.id, asset_digest=asset.digest, mode="archive", pushed=True)
    r.remote.link(ref=remote.id, asset_digest=asset.digest, mode="file", pushed=False)
    dest = tmp_path / "publish.sqlite"
    summary = export_publish_catalog(r, dest, HTTPS_PREFIX, **CATALOG_STORE)
    assert (summary["genome"], summary["asset"]) == (1, 1)
    return dest, genome_digest, asset


class TestExport:
    def test_export_contains_published_only(self, pushed):
        dest, genome_digest, asset = pushed
        eng = artifact_engine(dest)
        with eng.connect() as c:
            genomes = c.execute(select(Genome.__table__)).all()
            assert [g.digest for g in genomes] == [genome_digest]

            # aliases materialized from the (store-backed) manager into SQL rows
            aliases = c.execute(select(Alias.__table__)).all()
            assert [(a.name, a.genome_digest) for a in aliases] == [("rCRSd", genome_digest)]

            # only the pushed (archive) pair survives; local path stripped
            assets = c.execute(select(Asset.__table__)).all()
            assert [a.digest for a in assets] == [asset.digest]
            assert assets[0].path is None
            staged = c.execute(select(StagedAsset.__table__)).all()
            assert [(s.asset_digest, s.mode) for s in staged] == [(asset.digest, "archive")]
            assert staged[0].download_count == 0

            # one synthesized https remote; pushed link re-targeted at it
            remotes = c.execute(select(Remote.__table__)).all()
            assert len(remotes) == 1
            assert remotes[0].type == RemoteType.https
            assert remotes[0].prefix == HTTPS_PREFIX
            links = c.execute(select(RemoteAssetLink.__table__)).all()
            assert [(x.remote_id, x.asset_digest, x.mode, x.pushed) for x in links] == [
                (1, asset.digest, "archive", True)
            ]

            stamp = c.execute(select(AlembicVersion.__table__.c.version_num)).scalar()
            assert stamp == TARGET_ALEMBIC_VERSION


class TestImport:
    @pytest.fixture
    def artifact(self, pushed):
        return pushed

    @pytest.fixture
    def server_engine(self):
        eng = make_engine()
        SQLModel.metadata.create_all(eng)
        return eng

    def test_import_populates_server_catalog(self, artifact, server_engine):
        dest, genome_digest, asset = artifact
        summary = import_publish_catalog(server_engine, str(dest))
        assert summary["genome"] == 1

        with server_engine.connect() as c:
            assert c.execute(select(Genome.__table__.c.digest)).scalars().all() == [genome_digest]
            assert c.execute(select(Alias.__table__.c.name)).scalars().all() == ["rCRSd"]

    def test_import_yields_https_download_url(self, artifact, server_engine):
        """The redirect query the server runs must resolve to the https mirror."""
        dest, genome_digest, asset = artifact
        import_publish_catalog(server_engine, str(dest))

        from sqlmodel import Session

        from refgenie.db.tables import Configuration
        from refgenie.managers import RemoteManager

        with Session(server_engine) as session:
            stage_folder = session.scalar(select(Configuration.genome_stage_folder))
        tarball = pathlib.Path(stage_folder) / genome_digest / "fasta" / f"{asset.digest}.tgz"
        [(_, _, url)] = RemoteManager(server_engine).pushed_urls(asset.digest, {"archive": tarball})
        assert url == f"{HTTPS_PREFIX}/{genome_digest}/fasta/{asset.digest}.tgz"

    def test_reimport_is_idempotent_and_preserves_download_count(self, artifact, server_engine):
        dest, genome_digest, asset = artifact
        import_publish_catalog(server_engine, str(dest))

        # server-side activity between imports
        with server_engine.begin() as c:
            c.execute(
                StagedAsset.__table__.update()
                .where(StagedAsset.__table__.c.asset_digest == asset.digest)
                .values(download_count=7)
            )

        import_publish_catalog(server_engine, str(dest))
        with server_engine.connect() as c:
            staged = c.execute(select(StagedAsset.__table__)).all()
            assert len(staged) == 1  # merged by natural key, not duplicated
            assert staged[0].download_count == 7  # server-owned counter kept
            assert c.execute(select(Genome.__table__)).all().__len__() == 1

    def test_import_refuses_schema_mismatch(self, artifact, server_engine, tmp_path):
        import shutil

        dest = shutil.copy(artifact[0], tmp_path / "mismatch.sqlite")
        eng = artifact_engine(dest)
        with eng.begin() as c:
            c.execute(AlembicVersion.__table__.update().values(version_num="bogus"))
        eng.dispose()
        with pytest.raises(ValueError, match="refusing to import"):
            import_publish_catalog(server_engine, str(dest))

    def test_imported_catalog_serves_genomes_endpoint(self, artifact, tmp_path):
        """End to end: a server app bound to the imported catalog lists the genome."""
        requires_server()

        eng = make_engine()
        rgc = Refgenie(database_engine=eng, suppress_migrations=True)
        rgc.database.init(genome_folder=tmp_path / "server_genomes")

        dest, genome_digest, asset = artifact
        import_publish_catalog(eng, str(dest))
        # The catalog's store URL is fake; with no enabled store the server
        # serves DB-backed seqcol and lists every genome.
        with eng.begin() as conn:
            conn.execute(Store.__table__.update().values(enabled=False))

        with make_server_client(rgc) as client:
            resp = client.get("/v4/genomes")
            assert resp.status_code == 200
            digests = [g["digest"] for g in resp.json()["results"]]
            assert digests == [genome_digest]


class TestDisabledStores:
    """A disabled store is not published, and neither is anything of its genomes."""

    @staticmethod
    def _node(tmp_path, fixtures_path):
        """A node with a built, staged, pushed asset on a genome in store 'plantref'."""
        r = Refgenie(database_engine=make_engine(), suppress_migrations=True)
        r.database.init(genome_folder=tmp_path / "genomes", genome_stage_folder=tmp_path / "stage")
        register_fasta(r, fixtures_path)
        r.genome.initialize_genome(
            fasta_file_path=fixtures_path / "rCRSd.fa", alias_names=["rCRSd"], description="x"
        )
        r.build.run(
            recipe_name="fasta", genome_alias="rCRSd", asset_group_name="fasta", asset_name="t"
        )
        digest = r.alias.resolve("rCRSd")
        asset = r.asset.get(genome_digest=digest, asset_group_name="fasta", asset_name="t")
        r.stage.create(
            asset=asset, genome_folder=r.genome_folder, genome_stage_folder=r.genome_stage_folder
        )
        remote = r.remote.add(type=RemoteType.s3, prefix="b/assets", name="s3", push_command="true")
        r.remote.link(ref=remote.id, asset_digest=asset.digest, mode="archive", pushed=True)
        r.store.add("plantref", "url://plantref", StoreType.remote, priority=50)
        r.genome.set_store_name(digest, "plantref")
        AliasManager(r.database_engine).federated_sync(
            "athal", digest, "plantref", {"plantref": 50}
        )
        r.store.add("vgp", "url://vgp", StoreType.remote, priority=20)
        r.genome.add(FEDERATED_DIGEST, "federated genome", [], store_name="vgp")
        return r, digest

    @staticmethod
    def _set_enabled(r, name, enabled):
        with r.database_engine.begin() as conn:
            conn.execute(
                Store.__table__.update()
                .where(Store.__table__.c.name == name)
                .values(enabled=enabled)
            )

    @staticmethod
    def _counts(path):
        from refgenie.db.tables import Alias, Asset, AssetGroup, SeekKey, StagedAsset

        with artifact_engine(path).connect() as conn:
            return {
                m.__tablename__: len(conn.execute(select(m.__table__)).all())
                for m in (Alias, Asset, AssetGroup, SeekKey, StagedAsset, RemoteAssetLink)
            }

    def test_disabled_store_is_omitted_with_its_genomes(self, tmp_path, fixtures_path, caplog):
        r, digest = self._node(tmp_path, fixtures_path)
        self._set_enabled(r, "plantref", False)
        dest = tmp_path / "publish.sqlite"

        with caplog.at_level("INFO"):
            summary = export_publish_catalog(r, dest, HTTPS_PREFIX, **CATALOG_STORE)

        assert set(_exported_stores(dest)) == {"vgp", "refgenie-main"}
        assert _exported_genome_stores(dest) == {FEDERATED_DIGEST: "vgp"}
        assert self._counts(dest) == {
            "alias": 0,
            "asset": 0,
            "assetgroup": 0,
            "seekkey": 0,
            "stagedasset": 0,
            "remoteassetlink": 0,
        }
        assert summary["genome"] == 1
        assert "disabled store 'plantref': omitting its 1 genome(s)" in caplog.text

        server = _new_server(tmp_path)
        import_publish_catalog(server.database_engine, str(dest))
        assert [s.name for s in server.store.enabled_stores()] == ["vgp", "refgenie-main"]
        assert not server.genome.exists(digest)

    def test_re_enabling_the_store_brings_them_back(self, tmp_path, fixtures_path):
        r, digest = self._node(tmp_path, fixtures_path)
        self._set_enabled(r, "plantref", False)
        export_publish_catalog(r, tmp_path / "off.sqlite", HTTPS_PREFIX, **CATALOG_STORE)
        self._set_enabled(r, "plantref", True)
        dest = tmp_path / "on.sqlite"

        export_publish_catalog(r, dest, HTTPS_PREFIX, **CATALOG_STORE)

        assert set(_exported_stores(dest)) == {"vgp", "plantref", "refgenie-main"}
        assert _exported_genome_stores(dest) == {digest: "plantref", FEDERATED_DIGEST: "vgp"}
        counts = self._counts(dest)
        assert counts["asset"] == counts["assetgroup"] == counts["stagedasset"] == 1
        assert counts["remoteassetlink"] == 1 and counts["seekkey"] > 0
        assert counts["alias"] >= 2  # rCRSd from the store, athal from SQL

    def test_a_store_missing_from_the_registry_still_fails(self, tmp_path, fixtures_path):
        r, _ = self._node(tmp_path, fixtures_path)
        with r.database_engine.begin() as conn:
            conn.execute(
                Genome.__table__.update()
                .where(Genome.__table__.c.digest == FEDERATED_DIGEST)
                .values(store_name="nowhere")
            )
        with pytest.raises(RuntimeError, match="'nowhere'"):
            export_publish_catalog(r, tmp_path / "p.sqlite", HTTPS_PREFIX, **CATALOG_STORE)
