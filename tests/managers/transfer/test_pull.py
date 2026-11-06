"""Pulling assets from a remote server: ``rgc.transfer.pull`` and the AssetPuller."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from refgenie.models import GenomeAlias
from refgenie.managers.transfer.puller import PullOptions
from tests.support.factories import (
    GENOME,
    MOCK_MIRROR_URL,
    MOCK_SERVER_URL,
    add_asset_from_files,
    fake_digest,
    make_engine,
    seed_species,
)
from tests.support.server import needs_server


class TestPullFromServer:
    """Full rgc.transfer.pull() exercised end-to-end against the REAL server app.

    ``server_client_world`` (tests/conftest.py) serves the server through the
    REAL app (create_app) over an in-process ASGI transport -- no hand-rolled
    routes.
    """

    pytestmark = [pytest.mark.component, needs_server]

    def test_pull_asset_from_server(self, server_client_world):
        """Pulling a FASTA asset from a real remote server app registers it locally."""
        client_rg, _, url = server_client_world
        result = client_rg.transfer.pull(
            genome=GenomeAlias("rCRSd"),
            asset_group_name="fasta",
            force=True,
            force_large=True,
            force_server_urls=[url],
        )

        assert result is not None, "pull() returned None - asset was not pulled"
        assert client_rg.asset.exists(
            genome_digest=client_rg.alias.resolve("rCRSd"),
            asset_group_name="fasta",
            asset_name="default",
        )
        assert client_rg.alias.exists("rCRSd")
        assert client_rg.genome.exists(client_rg.alias.resolve("rCRSd"))

    def test_pulled_genome_keeps_the_servers_description(self, server_client_world):
        """A fasta pull keeps the server's genome description, even when forced.

        A forced re-pull removes the existing fasta and its genome, so the fasta
        step recreates the genome. It must take the description from what the
        server actually sends; reading a ``genome_description`` key the server
        never sends describes the recreated genome as "None".
        """
        client_rg, server_rg, url = server_client_world
        digest = server_rg.alias.resolve("rCRSd")
        server_description = server_rg.genome.get(digest).description
        assert server_description, "the server fixture must describe its genome"
        kwargs = dict(
            genome=GenomeAlias("rCRSd"),
            asset_group_name="fasta",
            force=True,
            force_large=True,
            force_server_urls=[url],
        )

        client_rg.transfer.pull(**kwargs)
        assert client_rg.genome.get(digest).description == server_description

        client_rg.transfer.pull(**kwargs)
        assert client_rg.genome.get(digest).description == server_description

    def test_pull_carries_the_servers_build_provenance(self, server_client_world):
        """A pulled name records the build the server recorded for it.

        Without this, ``build_digest`` is only addressable on the machine that
        built the asset, which defeats making builds addressable at all.
        """
        from tests.support.snapshots import asset_name_rows

        client_rg, server_rg, url = server_client_world
        served = {
            row.name: row for row in asset_name_rows(server_rg) if row.build_digest is not None
        }
        assert served, "the server fixture must have built something to serve"

        client_rg.transfer.pull(
            genome=GenomeAlias("rCRSd"),
            asset_group_name="fasta",
            force=True,
            force_large=True,
            force_server_urls=[url],
        )

        pulled = {row.name: row for row in asset_name_rows(client_rg)}
        for name, server_row in served.items():
            assert name in pulled, f"the server's name {name!r} was not adopted"
            assert pulled[name].build_digest == server_row.build_digest
            assert pulled[name].refgenie_version == server_row.refgenie_version
            assert pulled[name].build_level1 == server_row.build_level1
            # Surrogate ids are per-database; the server's numbering is not the
            # client's, so the recipe link is not carried over.
            assert pulled[name].recipe_id is None

    def test_pull_stops_after_first_server_succeeds(self, server_client_world):
        """With several servers, pull stops at the first that works, rather than
        pulling again from the next and failing with AssetExistsError."""
        client_rg, _, _ = server_client_world
        # force stays at its default (None), the case where a second pull raises.
        result = client_rg.transfer.pull(
            genome=GenomeAlias("rCRSd"),
            asset_group_name="fasta",
            force_large=True,
            force_server_urls=[MOCK_SERVER_URL, MOCK_MIRROR_URL],
        )

        assert result is not None
        assert client_rg.asset.exists(
            genome_digest=client_rg.alias.resolve("rCRSd"),
            asset_group_name="fasta",
            asset_name="default",
        )

    def test_pull_multiple_counts_multi_server_pull_as_success(self, server_client_world):
        """pull_multiple books a success when several servers are subscribed."""
        client_rg, server_rg, _ = server_client_world
        urls = [MOCK_SERVER_URL, MOCK_MIRROR_URL]
        genome_digest = server_rg.alias.resolve("rCRSd")

        puller = client_rg.transfer._puller
        with patch.object(client_rg.servers, "subscriptions", return_value=urls):
            successful = puller.pull_multiple(
                asset_list=[
                    {
                        "genome_digest": genome_digest,
                        "asset_group_name": "fasta",
                        "asset_name": "default",
                    }
                ],
                options=PullOptions(force_large=True),
            )

        assert len(successful.pulled) == 1, "multi-server pull should be counted as a success"

    def test_existing_asset_raises_without_force(self, server_client_world):
        """Pulling an already-present asset raises without force, replaces with it."""
        from refgenie.exceptions import AssetExistsError

        client_rg, _, url = server_client_world
        kwargs = dict(
            genome=GenomeAlias("rCRSd"),
            asset_group_name="fasta",
            force_large=True,
            force_server_urls=[url],
        )

        assert client_rg.transfer.pull(force=True, **kwargs) is not None

        with pytest.raises(AssetExistsError):
            client_rg.transfer.pull(force=False, **kwargs)

        # force=True over the existing asset re-downloads and replaces it.
        assert client_rg.transfer.pull(force=True, **kwargs) is not None
        assert client_rg.asset.exists(
            genome_digest=client_rg.alias.resolve("rCRSd"),
            asset_group_name="fasta",
            asset_name="default",
        )

    def test_declined_large_archive_is_skipped_not_failed(self, server_client_world, caplog):
        """Declining an archive over the size cutoff is a skip (None), not a RuntimeError,
        and the warning names the flags that override it, not --force."""
        client_rg, _, url = server_client_world
        result = client_rg.transfer.pull(
            genome=GenomeAlias("rCRSd"),
            asset_group_name="fasta",
            force_large=None,
            size_cutoff=0,
            confirm=lambda _message: False,
            force_server_urls=[url],
        )
        assert result is None
        assert not client_rg.alias.exists("rCRSd")
        assert "--pull-large" in caplog.text
        assert "--force" not in caplog.text

    def test_pull_multiple_counts_a_skip_as_skipped(self, server_client_world, caplog):
        """A skipped asset is booked as skipped in the summary, not as failed."""
        client_rg, server_rg, url = server_client_world
        entry = {
            "genome_digest": server_rg.alias.resolve("rCRSd"),
            "asset_group_name": "fasta",
            "asset_name": "default",
            "server_url": url,
        }
        result = client_rg.transfer._puller.pull_multiple(
            [entry], PullOptions(force_large=False, size_cutoff=0)
        )
        assert result.skipped == [entry] and not result.pulled
        assert "0 succeeded, 1 skipped, 0 failed" in caplog.text

    def test_a_genome_the_server_cannot_describe_is_named_as_a_server_fault(
        self, server_client_world
    ):
        """The server lists the alias but cannot return the genome's sequence
        collection, as api.refgenie.org did for t7 while it did not load the
        store holding it. The pull must say so, not that the alias is unknown."""
        from refgenie.exceptions import PullFailedError

        client_rg, server_rg, url = server_client_world
        digest = fake_digest("in_no_store")
        server_rg.genome.add(digest=digest, description="In no store", alias_names=["t7"])

        with pytest.raises(PullFailedError) as info:
            client_rg.transfer.pull(
                genome=GenomeAlias("t7"), asset_group_name="fasta", force_server_urls=[url]
            )
        reason = info.value.reasons[url]
        assert f"lists genome 't7' ({digest})" in reason
        assert "did not return its sequence collection" in reason
        assert "Could not resolve genome alias" not in reason
        assert not client_rg.alias.exists("t7")

    def test_an_alias_the_server_does_not_list_is_reported_unknown(self, server_client_world):
        from refgenie.exceptions import PullFailedError

        client_rg, _, url = server_client_world
        with pytest.raises(PullFailedError) as info:
            client_rg.transfer.pull(
                genome=GenomeAlias("no_such_genome"),
                asset_group_name="fasta",
                force_server_urls=[url],
            )
        assert info.value.reasons[url] == (
            f"Server {url} does not know genome alias 'no_such_genome'."
        )

    def test_pull_child_asset_records_its_parent(self, server_client_world, fixtures_path):
        """Pulling an asset that has a parent succeeds and links the parent locally,
        when both assets carry seek keys."""
        client_rg, server_rg, url = server_client_world
        genome_digest = server_rg.alias.resolve(GENOME)
        server_rg.asset_class.add(fixtures_path / "nested_seek_key_asset_class.yaml")
        client_rg.asset_class.add(fixtures_path / "nested_seek_key_asset_class.yaml")

        parent = server_rg.asset.get(
            genome_digest=genome_digest, asset_group_name="fasta", asset_name="default"
        )
        child = add_asset_from_files(
            server_rg,
            rel_dir="child_staging",
            files=["sub/inner/index.txt"],
            asset_class_name="nested_seek_key",
            asset_group_name="nested_seek_key",
            asset_name="default",
        )
        for group in ("fasta", "nested_seek_key"):
            assert server_rg.asset.seek_key.list_all(genome_digest, group, "default"), (
                f"the {group} asset must carry seek keys"
            )
        server_rg.asset.links.set_parents(
            genome_digest=genome_digest,
            asset_group_name="nested_seek_key",
            asset_name="default",
            parent_asset_digests=[parent.digest],
        )
        server_rg.stage.create(child, server_rg.genome_folder, server_rg.genome_stage_folder)

        common = dict(
            genome=GenomeAlias(GENOME), force=True, force_large=True, force_server_urls=[url]
        )
        assert client_rg.transfer.pull(asset_group_name="fasta", **common) is not None
        assert client_rg.transfer.pull(asset_group_name="nested_seek_key", **common) is not None

        parents = client_rg.asset.links.parents(
            genome_digest=genome_digest, asset_group_name="nested_seek_key", asset_name="default"
        )
        assert [p.digest for p in parents] == [parent.digest]

    def test_pull_carries_server_non_path_seek_key_values(self, server_client_world, fixtures_path):
        """A pulled asset keeps the server's string, json and extra seek key values,
        rather than re-deriving them from the local asset class."""
        import json

        client_rg, server_rg, url = server_client_world
        genome_digest = server_rg.alias.resolve(GENOME)
        server_rg.asset_class.add(fixtures_path / "test_asset_class.yaml")
        client_rg.asset_class.add(fixtures_path / "test_asset_class.yaml")
        server_values = {
            "test_version": "1.2.3",
            "test_metadata": json.dumps({"tool": "bwa", "threads": 4}),
            "version": "0.7.17",
        }
        asset = add_asset_from_files(
            server_rg,
            rel_dir="test_staging",
            files=["{genome}.extension", "{genome}.extension1", "{genome}.extension2"],
            asset_class_name="test",
            asset_group_name="test",
            asset_name="default",
            custom_seek_keys=server_values,
        )
        server_rg.stage.create(asset, server_rg.genome_folder, server_rg.genome_stage_folder)

        assert (
            client_rg.transfer.pull(
                genome=GenomeAlias(GENOME),
                asset_group_name="test",
                force=True,
                force_large=True,
                force_server_urls=[url],
            )
            is not None
        )

        for name, value in server_values.items():
            assert client_rg.asset.seek(genome_digest, "test", "default", seek_key_name=name) == (
                value
            ), f"pulled seek key {name!r} must carry the server's value"
        local = client_rg.asset.seek(genome_digest, "test", "default", seek_key_name="extension")
        assert Path(local).is_file()


#: An FHR record for the server's genome (only the fields refgenie promotes).
_WORM_FHR = {
    "schema": "https://raw.githubusercontent.com/FAIR-bioHeaders/FHR-Specification/main/fhr.json",
    "schemaVersion": 1,
    "genome": GENOME,
    "taxon": {"name": "Caenorhabditis elegans", "uri": "https://identifiers.org/taxonomy:6239"},
    "version": "3",
    "genomeSynonym": ["WBcel235"],
    "assemblySource": "Ensembl",
    "accessionID": {"name": "GCA_000002985.3", "url": "https://www.ncbi.nlm.nih.gov/"},
    "documentation": "worm assembly",
}


def _stage(server_rg, asset):
    server_rg.stage.create(asset, server_rg.genome_folder, server_rg.genome_stage_folder)


class TestPullOnAFreshInstall:
    """A client with no asset classes and no parents pulls straight from a server.

    The class definition comes from the server being pulled from, and missing
    parents are pulled first, from the same server. Real server app, no HTTP
    mocking.
    """

    pytestmark = [pytest.mark.component, needs_server]

    def test_the_asset_class_is_registered_from_the_server(self, bare_client_world):
        client_rg, server_rg, url, parent, _ = bare_client_world
        _stage(server_rg, parent)

        pulled = client_rg.transfer.pull(
            genome=GenomeAlias(GENOME), asset_group_name="fasta", force_server_urls=[url]
        )

        assert pulled is not None and pulled.digest == parent.digest
        local, remote = client_rg.asset_class.get("fasta"), server_rg.asset_class.get("fasta")
        assert local.version == remote.version
        assert local.serving_modes == remote.serving_modes
        assert {(k.name, k.value, k.type) for k in local.seek_keys} == {
            (k.name, k.value, k.type) for k in remote.seek_keys
        }, "the class must arrive with its seek keys, so it can build locally too"
        assert Path(client_rg.asset.seek(parent.genome_digest, "fasta", "default")).is_file()

    def test_a_missing_parent_is_pulled_first(self, bare_client_world):
        client_rg, server_rg, url, parent, child = bare_client_world
        _stage(server_rg, parent)
        _stage(server_rg, child)
        genome_digest = server_rg.alias.resolve(GENOME)

        pulled = client_rg.transfer.pull(
            genome=GenomeAlias(GENOME),
            asset_group_name="child_index",
            asset_name="v1",
            force_server_urls=[url],
        )

        assert pulled is not None and pulled.digest == child.digest
        assert client_rg.asset.exists_digest(parent.digest), "parent not pulled"
        parents = client_rg.asset.links.parents(
            genome_digest=genome_digest, asset_group_name="child_index", asset_name="v1"
        )
        assert [p.digest for p in parents] == [parent.digest]
        assert client_rg.asset_class.get("nested_seek_key").name == "nested_seek_key"
        assert client_rg.asset_class.get("fasta").name == "fasta"

    def test_a_parent_that_cannot_be_pulled_is_named(self, bare_client_world):
        from refgenie.exceptions import ParentPullError

        client_rg, server_rg, url, parent, child = bare_client_world
        _stage(server_rg, child)  # the parent is never staged: it cannot be served
        genome_digest = server_rg.alias.resolve(GENOME)

        with pytest.raises(ParentPullError) as excinfo:
            client_rg.transfer.pull(
                genome=GenomeAlias(GENOME),
                asset_group_name="child_index",
                asset_name="v1",
                force_server_urls=[url],
            )

        message = str(excinfo.value)
        for part in (genome_digest, "fasta", "default", parent.digest, "child_index:v1"):
            assert part in message, f"{part!r} missing from: {message}"
        assert "No archive found" in message, f"the reason is missing from: {message}"
        # Rollback removed only what this pull created: the genome and alias it
        # registered, and nothing of the child.
        assert not client_rg.asset.exists_digest(child.digest)
        assert not client_rg.genome.exists(genome_digest)
        assert not client_rg.alias.exists(GENOME)

    def test_a_pulled_parent_survives_a_failed_child(self, bare_client_world):
        from refgenie.exceptions import PullFailedError

        client_rg, server_rg, url, parent, child = bare_client_world
        _stage(server_rg, parent)  # the child is never staged: it cannot be served
        genome_digest = server_rg.alias.resolve(GENOME)

        with pytest.raises(PullFailedError):
            client_rg.transfer.pull(
                genome=GenomeAlias(GENOME),
                asset_group_name="child_index",
                asset_name="v1",
                force_server_urls=[url],
            )

        # The parent is a complete asset of its own, so the child's rollback
        # must not take it, or the genome it lives in, away.
        assert client_rg.asset.exists_digest(parent.digest)
        assert client_rg.genome.exists(genome_digest)
        assert client_rg.alias.resolve(GENOME) == genome_digest
        assert not client_rg.asset.exists_digest(child.digest)

    def test_the_servers_genome_record_is_copied(self, bare_client_world):
        """A pull that creates the genome copies the server's species, taxon
        and assembly fields, and its FHR record, not just the description."""
        client_rg, server_rg, url, parent, child = bare_client_world
        _stage(server_rg, parent)
        _stage(server_rg, child)
        digest = server_rg.alias.resolve(GENOME)
        server_rg.genome.apply_fhr(digest, _WORM_FHR)
        seed_species(server_rg, 6239, "Caenorhabditis elegans", "worm", ["nematode worm"])

        # A non-fasta asset first: the genome is created by the pull's genome
        # resolution and the fasta arrives as an auto-pulled parent.
        client_rg.transfer.pull(
            genome=GenomeAlias(GENOME),
            asset_group_name="child_index",
            asset_name="v1",
            force_server_urls=[url],
        )

        local = client_rg.genome.get(digest)
        assert local.species_name == "Caenorhabditis elegans"
        assert local.synonyms == ["WBcel235"]
        assert local.taxon_id == 6239
        assert local.assembly_source == "Ensembl"
        assert local.assembly_accession == "GCA_000002985.3"
        assert local.description == "worm assembly"
        fhr = client_rg.refget_store.get_fhr_metadata(digest)
        assert fhr is not None and fhr.to_dict()["genome"] == GENOME
        species = client_rg.species.get(6239)
        assert (species.scientific_name, species.common_name) == ("Caenorhabditis elegans", "worm")
        assert species.synonyms == ["nematode worm"]
        assert client_rg.genome.detail(digest).common_name == "worm"

    def test_a_later_pull_replaces_the_species_row(self, bare_client_world):
        client_rg, server_rg, url, parent, child = bare_client_world
        _stage(server_rg, parent)
        _stage(server_rg, child)
        digest = server_rg.alias.resolve(GENOME)
        server_rg.genome.apply_fhr_columns(digest, {"taxon": {"uri": "taxonomy:6239"}})
        seed_species(server_rg, 6239, "Caenorhabditis elegans", "roundworm", ["nematode"])
        client_rg.transfer.pull(
            genome=GenomeAlias(GENOME), asset_group_name="fasta", force_server_urls=[url]
        )
        assert client_rg.species.get(6239).synonyms == ["nematode"]

        seed_species(server_rg, 6239, "Caenorhabditis elegans", "C. elegans", ["worm"])
        client_rg.transfer.pull(
            genome=GenomeAlias(GENOME),
            asset_group_name="child_index",
            asset_name="v1",
            force_server_urls=[url],
        )

        species = client_rg.species.get(6239)
        assert (species.common_name, species.synonyms) == ("C. elegans", ["worm"])

    def test_local_genome_values_are_not_overwritten(self, bare_client_world):
        """Fields the user set locally survive a pull; unset ones are filled."""
        client_rg, server_rg, url, parent, child = bare_client_world
        _stage(server_rg, parent)
        _stage(server_rg, child)
        digest = server_rg.alias.resolve(GENOME)
        server_rg.genome.apply_fhr_columns(digest, {"taxon": {"name": "Caenorhabditis elegans"}})

        client_rg.transfer.pull(
            genome=GenomeAlias(GENOME), asset_group_name="fasta", force_server_urls=[url]
        )
        assert client_rg.genome.get(digest).species_name == "Caenorhabditis elegans"

        client_rg.genome.apply_fhr_columns(digest, {"taxon": {"name": "my worm"}})
        server_rg.genome.apply_fhr_columns(digest, _WORM_FHR)
        seed_species(server_rg, 6239, "Caenorhabditis elegans", "roundworm")
        client_rg.transfer.pull(
            genome=GenomeAlias(GENOME),
            asset_group_name="child_index",
            asset_name="v1",
            force_server_urls=[url],
        )

        local = client_rg.genome.get(digest)
        assert local.species_name == "my worm", "a local value was overwritten"
        assert local.taxon_id == 6239, "an unset local value was not filled"
        assert client_rg.species.get(6239).common_name == "roundworm"
        assert client_rg.genome.detail(digest).common_name == "roundworm"


# ===========================================================================
# File-mode download + serving-mode dispatch
# ===========================================================================


_MOCK_ALIASES = [{"name": "test_genome", "genome_digest": "genome_digest_001"}]


class TestPullFileMode:
    """Tests for the _pull_file_mode file-by-file download helper."""

    pytestmark = pytest.mark.component

    def test_downloads_files_into_nested_dirs(self, refgenie_fs, tmp_path):
        """_pull_file_mode downloads each file from the list, creating subdirs."""
        from refgenie.managers.transfer.puller import PullTransaction

        puller = refgenie_fs.transfer._puller
        client = MagicMock()
        client.asset.files.return_value = ["genome.gtf", "subdir/annotation.bed"]

        def mock_download(asset_digest, file_path, output_path):
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text(f"content of {file_path}")
            return output_path

        client.asset.download_file.side_effect = mock_download

        asset_dir = tmp_path / "test_asset"
        puller._pull_file_mode(
            client=client,
            asset_digest="abc123",
            asset_dir=asset_dir,
            bundle_name="test/fasta:default",
            txn=PullTransaction(puller._alias_manager, puller._genome_manager),
        )

        assert client.asset.download_file.call_count == 2
        assert (asset_dir / "genome.gtf").exists()
        assert (asset_dir / "subdir" / "annotation.bed").exists()

    def test_empty_file_list_raises(self, refgenie_fs, tmp_path):
        """_pull_file_mode raises RuntimeError (never downloads) on an empty list."""
        from refgenie.managers.transfer.puller import PullTransaction

        puller = refgenie_fs.transfer._puller
        client = MagicMock()
        client.asset.files.return_value = []

        with pytest.raises(RuntimeError, match="No files found"):
            puller._pull_file_mode(
                client=client,
                asset_digest="abc123",
                asset_dir=tmp_path / "test_asset",
                bundle_name="test/fasta:default",
                txn=PullTransaction(puller._alias_manager, puller._genome_manager),
            )
        client.asset.download_file.assert_not_called()

    def test_rollback_removes_created_dir_on_failure(self, refgenie_fs, tmp_path):
        """On a mid-download failure, PullTransaction removes the created directory."""
        from refgenie.managers.transfer.puller import PullTransaction

        puller = refgenie_fs.transfer._puller
        client = MagicMock()
        client.asset.files.return_value = ["file1.txt", "file2.txt"]

        call_count = [0]

        def mock_download(asset_digest, file_path, output_path):
            call_count[0] += 1
            if call_count[0] == 1:
                output_path.parent.mkdir(parents=True, exist_ok=True)
                output_path.write_text("content")
                return output_path
            raise ConnectionError("Download failed")

        client.asset.download_file.side_effect = mock_download

        asset_dir = tmp_path / "rollback_test"
        txn = PullTransaction(puller._alias_manager, puller._genome_manager)
        with pytest.raises(ConnectionError):
            with txn:
                puller._pull_file_mode(
                    client=client,
                    asset_digest="abc123",
                    asset_dir=asset_dir,
                    bundle_name="test/fasta:default",
                    txn=txn,
                )

        assert asset_dir in txn._created_dirs  # tracked for rollback
        assert not asset_dir.exists()  # and actually removed


class TestFinalizePulledAsset:
    pytestmark = pytest.mark.component

    def test_nested_asset_dirs_appear_under_the_alias_tree(self, refgenie_fs, fixtures_path):
        from refgenie.managers.transfer.puller import PulledAssetMetadata

        r = refgenie_fs
        r.asset_class.add(fixtures_path / "nested_seek_key_asset_class.yaml")
        digest = r.alias.resolve(GENOME)

        asset_dir = r.genome_folder / "pulled_nested_asset"
        (asset_dir / "sub" / "inner").mkdir(parents=True)

        meta = PulledAssetMetadata(
            asset_digest="f" * 64,
            declared_modes=["file"],
            asset_parents=[],
            asset_name="pulled",
            asset_class_name="nested_seek_key",
            asset_names=[{"name": "pulled"}],
        )

        added = r.transfer._puller._finalize_pulled_asset(
            alias_name=GENOME,
            genome_digest=digest,
            asset_group_name="nested_group",
            asset_name="pulled",
            asset_dir=asset_dir,
            meta=meta,
            genome_attrs={},
        )

        assert added is not None
        alias_dir = r.alias_folder / GENOME / "nested_group" / "pulled"
        assert (alias_dir / "sub" / "inner").is_dir()
        assert r.asset.group.get_default("nested_group", genome_digest=digest) == "pulled"


def _serve_with_modes(engine, tmp_path, fixtures_path, modes):
    """A served world whose fasta asset class offers ``modes``, staged accordingly."""
    from sqlmodel import Session, select

    from refgenie.db.tables import AssetClass
    from tests.support.server import make_server_client_world, serve_refgenie

    server_rg, client_rg = make_server_client_world(engine, tmp_path, fixtures_path, stage=False)
    with Session(server_rg.database_engine) as session:
        asset_class = session.exec(select(AssetClass).where(AssetClass.name == "fasta")).first()
        asset_class.serving_modes = modes
        session.add(asset_class)
        session.commit()
    asset = next(iter(server_rg.asset.list_assets()))
    server_rg.stage.create(asset, server_rg.genome_folder, server_rg.genome_stage_folder)
    return client_rg, serve_refgenie(client_rg, server_rg, MOCK_SERVER_URL)


class TestPullDownloadModeSelect:
    """pull() downloads the way the server serves the asset."""

    pytestmark = [pytest.mark.component, needs_server]

    @pytest.mark.parametrize(
        "modes, expected",
        [(["archive", "file"], "archive"), (["file"], "file")],
        ids=["both_prefer_archive", "only_file"],
    )
    def test_selects_download_mode(self, engine, tmp_path, fixtures_path, modes, expected):
        client_rg, served = _serve_with_modes(engine, tmp_path, fixtures_path, modes)
        puller = client_rg.transfer._puller
        with (
            served,
            patch.object(puller, "_pull_archive_mode", wraps=puller._pull_archive_mode) as arc,
            patch.object(puller, "_pull_file_mode", wraps=puller._pull_file_mode) as file,
        ):
            client_rg.transfer.pull(
                genome=GenomeAlias(GENOME),
                asset_group_name="fasta",
                force_large=True,
                force_server_urls=[MOCK_SERVER_URL],
            )
        assert (arc.called, file.called) == (expected == "archive", expected == "file")
        assert client_rg.asset.exists(
            genome_digest=client_rg.alias.resolve(GENOME),
            asset_group_name="fasta",
            asset_name="default",
        )

    def test_none_mode_raises_pull_failed(self, engine, tmp_path, fixtures_path):
        """A metadata-only asset cannot be pulled."""
        from refgenie.exceptions import PullFailedError

        client_rg, served = _serve_with_modes(engine, tmp_path, fixtures_path, ["none"])
        with served, pytest.raises(PullFailedError, match="metadata-only"):
            client_rg.transfer.pull(
                genome=GenomeAlias(GENOME),
                asset_group_name="fasta",
                force_server_urls=[MOCK_SERVER_URL],
            )


# ===========================================================================
# Failure surfacing (mock-only, unit tier)
# ===========================================================================


@pytest.fixture
def mock_puller(tmp_path):
    """A real AssetPuller over MagicMock collaborators, with one subscription."""
    from refgenie.exceptions import MissingAssetGroupError
    from refgenie.managers.transfer.puller import AssetPuller

    engine = make_engine()
    servers = MagicMock()
    asset = MagicMock()
    alias_manager = MagicMock()
    genome_manager = MagicMock()
    links = MagicMock()
    groups = MagicMock()

    asset.genome_folder = tmp_path / "genomes"
    asset.data_folder = tmp_path / "genomes" / "data"
    asset.alias_folder = tmp_path / "aliases"
    asset.links = links
    asset.group = groups

    puller = AssetPuller(
        database_engine=engine,
        asset=asset,
        servers=servers,
        alias_manager=alias_manager,
        genome_manager=genome_manager,
        asset_class=MagicMock(),
    )
    servers.subscriptions.return_value = ["http://test-server:5000"]
    alias_manager.resolve.return_value = "abc123digest"
    alias_manager.exists.return_value = True
    # A MagicMock is truthy, which would make every token a known digest.
    genome_manager.exists.return_value = False
    servers.client.return_value.alias.list.return_value = []
    asset.exists.return_value = False
    groups.get_default.side_effect = MissingAssetGroupError("genome", "fasta")
    return puller


class TestResolvePullGenome:
    """``pull`` takes a ``GenomeAlias`` or a ``GenomeDigest``; its type says
    which, and the puller works out which alias the pull is addressed to."""

    def _resolve(self, puller, genome):
        from refgenie.managers.transfer.puller import PullTransaction

        client = puller._servers.client("http://test-server:5000")
        with PullTransaction(puller._alias_manager, puller._genome_manager) as txn:
            # The third item is the server's genome record; these tests
            # check only the alias and digest the pull is addressed to.
            return puller._resolve_pull_genome(genome, client, "http://test-server:5000", txn)[:2]

    def test_a_local_alias_is_used_as_given(self, mock_puller):
        assert self._resolve(mock_puller, GenomeAlias("hg38")) == ("hg38", "abc123digest")
        mock_puller._genome_manager.ensure_from_remote.assert_not_called()

    def test_a_digest_is_addressed_to_the_servers_first_alias(self, mock_puller):
        digest = fake_digest("d")
        mock_puller._genome_manager.exists.return_value = True
        mock_puller._servers.client.return_value.alias.list.return_value = [
            {"name": "first", "genome_digest": digest},
            {"name": "second", "genome_digest": digest},
        ]
        assert self._resolve(mock_puller, digest) == ("first", digest)
        assert mock_puller._genome_manager.ensure_from_remote.call_count == 2
        mock_puller._alias_manager.resolve.assert_not_called()

    def test_an_unknown_alias_is_registered_from_the_server(self, mock_puller):
        from refgenie.exceptions import MissingAliasError
        from refgenie.managers.genome import GenomeCreationResult

        digest = fake_digest("d")
        mock_puller._alias_manager.resolve.side_effect = [MissingAliasError("remote"), digest]
        ensure = mock_puller._genome_manager.ensure_from_remote
        ensure.return_value = GenomeCreationResult(
            success=True, genome_digest=digest, alias_name=GenomeAlias("remote")
        )
        assert self._resolve(mock_puller, GenomeAlias("remote")) == ("remote", digest)
        assert ensure.call_args.kwargs["alias_name"] == "remote"
        # An alias is never tried as a digest.
        mock_puller._servers.client.return_value.alias.list.assert_not_called()

    def test_an_alias_the_server_does_not_know_cannot_be_served(self, mock_puller):
        from refgenie.exceptions import MissingAliasError, ServerCannotServe
        from refgenie.managers.genome import GenomeCreationResult

        mock_puller._alias_manager.resolve.side_effect = MissingAliasError("nope")
        mock_puller._genome_manager.ensure_from_remote.return_value = GenomeCreationResult(
            success=False
        )
        with pytest.raises(ServerCannotServe, match="genome alias 'nope'"):
            self._resolve(mock_puller, GenomeAlias("nope"))

    def test_a_digest_with_no_server_alias_falls_back_to_a_local_one(self, mock_puller):
        digest = fake_digest("d")
        mock_puller._genome_manager.exists.return_value = True
        mock_puller._alias_manager.get_for_genome.return_value = ["local"]
        assert self._resolve(mock_puller, digest) == ("local", digest)

    def test_a_plain_string_is_refused(self, mock_puller):
        """A plain str could be either kind; pull makes the caller say which."""
        with pytest.raises(TypeError, match="GenomeAlias or a GenomeDigest"):
            mock_puller.pull(genome="hg38", asset_group_name="fasta")
        mock_puller._servers.client.assert_not_called()


class TestPullFailures:
    """pull() raises, naming each server's reason, when servers answer badly."""

    @pytest.mark.parametrize(
        "exc",
        [ValueError("Missing operation ID"), ConnectionError("Connection refused")],
        ids=["value_error", "connection_error"],
    )
    def test_unreachable_asset_groups_raises_with_reason(self, mock_puller, caplog, exc):
        """When asset_group.list raises, the single server is exhausted; pull logs
        the failure and raises PullFailedError naming that server's reason."""
        from refgenie.exceptions import NoSubscriptionsError, PullFailedError

        client = MagicMock()
        client.asset_group.list.side_effect = exc
        mock_puller._servers.client.return_value = client

        with pytest.raises(PullFailedError) as info:
            mock_puller.pull(genome=GenomeAlias("testgenome"), asset_group_name="fasta")
        assert not isinstance(info.value, NoSubscriptionsError)
        reason = info.value.reasons["http://test-server:5000"]
        assert "Failed to query asset groups from http://test-server:5000" in reason
        assert str(exc) in str(info.value)
        assert "Failed to query asset groups from http://test-server:5000" in caplog.text

    def test_unreachable_relationships_raises_with_reason(self, mock_puller, caplog):
        """When the relationships query raises (a 500 from the server, say), the
        server is skipped like any other unreachable query, and the failure the
        caller sees says so rather than blaming the subscriptions."""
        from refgenie.exceptions import NoSubscriptionsError, PullFailedError

        client = MagicMock()
        client.asset_group.list.return_value = [
            {"id": "group1", "name": "fasta", "genome_digest": "abc123digest"}
        ]
        client.asset.list.return_value = [
            {"id": "asset1", "name": "default", "digest": "d1", "serving_modes": []}
        ]
        client.asset.relationships.side_effect = RuntimeError(
            "500 Server Error: Internal Server Error"
        )
        mock_puller._servers.client.return_value = client

        with pytest.raises(PullFailedError) as info:
            mock_puller.pull(genome=GenomeAlias("testgenome"), asset_group_name="fasta")
        assert not isinstance(info.value, NoSubscriptionsError)
        assert list(info.value.reasons) == ["http://test-server:5000"]
        assert "Failed to query asset relationships" in str(info.value)
        assert "500 Server Error" in str(info.value)
        assert "Failed to query asset relationships from http://test-server:5000" in caplog.text

    def test_every_server_reason_is_collected(self, mock_puller):
        """With two servers failing, both reasons are carried, in order."""
        from refgenie.exceptions import PullFailedError

        mock_puller._servers.subscriptions.return_value = ["http://a", "http://b"]
        client = MagicMock()
        client.asset_group.list.side_effect = ConnectionError("Connection refused")
        mock_puller._servers.client.return_value = client

        with pytest.raises(PullFailedError) as info:
            mock_puller.pull(genome=GenomeAlias("testgenome"), asset_group_name="fasta")
        assert list(info.value.reasons) == ["http://a", "http://b"]
        assert "http://a" in str(info.value) and "http://b" in str(info.value)

    def test_zero_asset_groups_raises(self, mock_puller):
        """Zero asset groups returned raises PullFailedError."""
        from refgenie.exceptions import PullFailedError

        client = MagicMock()
        client.asset_group.list.return_value = []
        mock_puller._servers.client.return_value = client

        with pytest.raises(PullFailedError, match="Expected one asset group, got 0"):
            mock_puller.pull(genome=GenomeAlias("testgenome"), asset_group_name="fasta")

    def test_multiple_assets_no_name_raises_helpful(self, mock_puller):
        """Multiple assets with no name specified raises with a 'specify one' hint."""
        from refgenie.exceptions import PullFailedError

        client = MagicMock()
        client.asset_group.list.return_value = [
            {"id": "group1", "name": "fasta", "genome_digest": "abc123digest"}
        ]
        client.asset.list.return_value = [
            {"id": "asset1", "name": "default", "digest": "d1"},
            {"id": "asset2", "name": "v2", "digest": "d2"},
        ]
        mock_puller._servers.client.return_value = client

        with pytest.raises(PullFailedError, match="Multiple assets found.*Specify one with"):
            mock_puller.pull(
                genome=GenomeAlias("testgenome"), asset_group_name="fasta", asset_name=None
            )

    def test_a_client_side_error_is_not_blamed_on_the_server(self, mock_puller):
        """A connection error during alias resolution is reported as itself,
        even when the server's alias listing would name the genome."""
        import httpx

        from refgenie.exceptions import MissingAliasError, PullFailedError
        from refgenie.managers.genome import GenomeCreationResult

        mock_puller._alias_manager.resolve.side_effect = MissingAliasError("testgenome")
        mock_puller._genome_manager.ensure_from_remote.return_value = GenomeCreationResult(
            success=False
        )
        client = mock_puller._servers.client.return_value
        client.alias.resolve.side_effect = httpx.ConnectError("Name or service not known")
        client.alias.list.return_value = [{"name": "testgenome", "genome_digest": "d"}]

        with pytest.raises(PullFailedError, match="Name or service not known") as info:
            mock_puller.pull(genome=GenomeAlias("testgenome"), asset_group_name="fasta")
        assert "maintainers" not in str(info.value)

    def test_no_servers_raises(self, mock_puller):
        """No server subscriptions raises NoSubscriptionsError, a PullFailedError."""
        from refgenie.exceptions import NoSubscriptionsError

        mock_puller._servers.subscriptions.return_value = []
        with pytest.raises(NoSubscriptionsError, match="No server subscriptions found"):
            mock_puller.pull(genome=GenomeAlias("testgenome"), asset_group_name="fasta")


class TestArchiveVerification:
    """Archive verification must use the tarball byte-digest, or say it did not run."""

    def test_warns_when_tarball_digest_absent(self, mock_puller, tmp_path, caplog):
        """A staged record without ``tarball_digest`` warns that it could not
        verify, rather than comparing against the asset's identity digest."""
        import logging

        from refgenie.managers.transfer.puller import PullTransaction

        client = MagicMock()

        def fake_download(operation_id, output_path, url_format_params, name):
            Path(output_path).parent.mkdir(parents=True, exist_ok=True)
            Path(output_path).write_bytes(b"tarball bytes")

        client.download_with_progress.side_effect = fake_download

        staged = {"mode": "archive", "digest": "asset-identity-digest", "size": 999}
        with caplog.at_level(logging.WARNING):
            with patch("refgenie.managers.transfer.puller.untar"):
                mock_puller._pull_archive_mode(
                    client=client,
                    staged_asset_metadata=staged,
                    asset_digest="asset-identity-digest",
                    genome_digest="g1",
                    asset_group_name="fasta",
                    asset_name="default",
                    asset_dir=tmp_path / "data" / "g1" / "fasta" / "asset-identity-digest",
                    bundle_name="g1/fasta:default",
                    txn=PullTransaction(mock_puller._alias_manager, mock_puller._genome_manager),
                    options=PullOptions(force_large=True),
                )

        assert "could not be verified" in caplog.text.lower()


class TestPullSigintHandling:
    """The SIGINT handler installed around the archive download is restored afterwards."""

    def _run(self, mock_puller, tmp_path, *, download_side_effect=None):
        from refgenie.managers.transfer.puller import PullTransaction

        client = MagicMock()
        if download_side_effect is not None:
            client.download_with_progress.side_effect = download_side_effect
        else:

            def fake_download(operation_id, output_path, url_format_params, name):
                Path(output_path).parent.mkdir(parents=True, exist_ok=True)
                Path(output_path).write_bytes(b"tarball bytes")

            client.download_with_progress.side_effect = fake_download

        staged = {"mode": "archive", "tarball_digest": None, "tarball_size": 999}
        with patch("refgenie.managers.transfer.puller.untar"):
            mock_puller._pull_archive_mode(
                client=client,
                staged_asset_metadata=staged,
                asset_digest="asset-digest",
                genome_digest="g1",
                asset_group_name="fasta",
                asset_name="default",
                asset_dir=tmp_path / "data" / "g1" / "fasta" / "asset-digest",
                bundle_name="g1/fasta:default",
                txn=PullTransaction(mock_puller._alias_manager, mock_puller._genome_manager),
                options=PullOptions(force_large=True),
            )

    def test_sigint_handler_restored_after_successful_download(self, mock_puller, tmp_path):
        import signal

        original_handler = signal.getsignal(signal.SIGINT)
        try:
            self._run(mock_puller, tmp_path)
            assert signal.getsignal(signal.SIGINT) is original_handler
        finally:
            signal.signal(signal.SIGINT, original_handler)

    def test_sigint_handler_restored_after_failed_download(self, mock_puller, tmp_path):
        import signal

        original_handler = signal.getsignal(signal.SIGINT)
        try:
            with pytest.raises(RuntimeError):
                self._run(
                    mock_puller,
                    tmp_path,
                    download_side_effect=RuntimeError("download failed"),
                )
            assert signal.getsignal(signal.SIGINT) is original_handler
        finally:
            signal.signal(signal.SIGINT, original_handler)


# ===========================================================================
# Atomicity / rollback, against a server that fails
# ===========================================================================


def _break_archives(server_rg, how):
    """Delete or corrupt every staged archive on the server."""
    archives = list(Path(server_rg.genome_stage_folder).rglob("*.tgz"))
    assert archives
    for archive in archives:
        archive.unlink() if how == "delete" else archive.write_bytes(b"not the archive")


def _pull_rcrsd(client_rg, url):
    client_rg.transfer.pull(
        genome=GenomeAlias(GENOME),
        asset_group_name="fasta",
        force_large=True,
        force_server_urls=[url],
    )


def _files_under(folder):
    return [p for p in Path(folder).rglob("*") if p.is_file()]


class TestPullRollback:
    """A failed pull leaves no genome, alias or files behind, and keeps what was there."""

    pytestmark = [pytest.mark.component, needs_server]

    @pytest.mark.parametrize(
        "how, error",
        [("delete", "404 Not Found"), ("corrupt", "checksum mismatch")],
        ids=["download_fails", "checksum_mismatch"],
    )
    def test_failed_download_rolls_back(self, server_client_world, how, error):
        client_rg, server_rg, url = server_client_world
        digest = server_rg.alias.resolve(GENOME)
        _break_archives(server_rg, how)

        with pytest.raises(Exception, match=error):
            _pull_rcrsd(client_rg, url)

        assert not client_rg.alias.exists(GENOME)
        assert not client_rg.genome.exists(digest)
        assert _files_under(Path(client_rg.genome_folder) / "data") == []

    def test_rollback_when_server_has_no_archive(self, bare_client_world):
        from refgenie.exceptions import NoArchiveError

        client_rg, server_rg, url, _, _ = bare_client_world
        digest = server_rg.alias.resolve(GENOME)

        with pytest.raises(NoArchiveError):
            _pull_rcrsd(client_rg, url)

        assert not client_rg.alias.exists(GENOME)
        assert not client_rg.genome.exists(digest)
        assert not (Path(client_rg.genome_folder) / "data" / digest).exists()

    def test_rollback_preserves_genome_known_under_other_alias(self, server_client_world):
        """Only the alias the failed pull added is removed, never the genome."""
        client_rg, server_rg, url = server_client_world
        digest = server_rg.alias.resolve(GENOME)
        client_rg.genome.add(digest=digest, description="Local", alias_names=["old_alias"])
        _break_archives(server_rg, "delete")

        with pytest.raises(Exception, match="404"):
            _pull_rcrsd(client_rg, url)

        assert client_rg.genome.exists(digest)
        assert client_rg.alias.exists("old_alias")
        assert not client_rg.alias.exists(GENOME)

    def test_no_rollback_when_genome_preexists(self, server_client_world):
        client_rg, server_rg, url = server_client_world
        digest = server_rg.alias.resolve(GENOME)
        client_rg.genome.add(digest=digest, description="Local", alias_names=[GENOME])
        _break_archives(server_rg, "delete")

        with pytest.raises(Exception, match="404"):
            _pull_rcrsd(client_rg, url)

        assert client_rg.alias.exists(GENOME)
        assert client_rg.genome.exists(digest)


# ===========================================================================
# Bulk pull helpers (mock-only, unit)
# ===========================================================================


class TestBuildPullParents:
    """build.run(pull_parents=True) delegates to pull() with the right arguments."""

    def test_build_run_pull_parents_passes_genome_as_alias(self, refgenie_minimal, monkeypatch):
        """``--pull-parents`` pulls the *fasta* group for the *genome* alias."""
        from refgenie.exceptions import MissingAliasError

        calls = {}

        def fake_pull(*args, **kwargs):
            calls["args"] = args
            calls["kwargs"] = kwargs
            return None

        monkeypatch.setattr(refgenie_minimal.transfer, "pull", fake_pull)

        with pytest.raises(MissingAliasError):
            refgenie_minimal.build.run(
                recipe_name="fasta",
                genome_alias="nonexistent_genome",
                asset_group_name="fasta",
                pull_parents=True,
            )

        assert calls["args"] == ()
        assert calls["kwargs"] == {
            "asset_group_name": "fasta",
            "genome": "nonexistent_genome",
        }


class TestPullMultiple:
    """pull_multiple continues past individual failures and skips invalid entries."""

    def test_continues_on_individual_failure(self, refgenie_minimal):
        asset_list = [
            {
                "genome_digest": fake_digest("abc123"),
                "asset_group_name": "fasta",
                "asset_name": "default",
                "server_url": "http://example.com",
            },
            {
                "genome_digest": fake_digest("def456"),
                "asset_group_name": "fasta",
                "asset_name": "default",
                "server_url": "http://example.com",
            },
        ]

        call_count = 0

        def mock_pull(*args, **kwargs):
            nonlocal call_count
            call_count += 1
            if call_count == 1:
                raise Exception("Simulated failure")
            return MagicMock()

        with patch.object(refgenie_minimal.transfer._puller, "pull", side_effect=mock_pull):
            results = refgenie_minimal.transfer._puller.pull_multiple(asset_list)
        assert len(results.pulled) == 1
        assert len(results.failed) == 1

    def test_skips_entries_missing_required_fields(self, refgenie_minimal):
        asset_list = [
            {"genome_digest": None, "asset_group_name": "fasta"},
            {"genome_digest": "abc123"},  # missing asset_group_name
        ]
        assert refgenie_minimal.transfer._puller.pull_multiple(asset_list).failed == asset_list


class TestOverwriteConfirmation:
    """The "replace the existing asset directory?" prompt goes through the
    caller's confirmer, not a bare `rich.prompt.Confirm.ask`.

    So `confirm=deny` guards that branch. Calling rich directly there would
    leave it unguarded, and a pull driven from a server request or a worker
    thread would block forever on a terminal nobody is watching.
    """

    @pytest.fixture
    def puller_at_the_prompt(self, mock_puller):
        """A puller positioned exactly on the "directory already exists" branch."""
        from refgenie.managers.transfer.puller import PulledAssetMetadata

        asset_digest = "existingassetdigest"
        mock_puller._fetch_asset_metadata = MagicMock(
            return_value=PulledAssetMetadata(
                asset_digest=asset_digest,
                declared_modes=["archive"],
                asset_parents=[],
                asset_name="default",
                asset_class_name="fasta",
                asset_names=[],
            )
        )
        mock_puller.get_client = MagicMock(return_value=MagicMock())
        asset_dir = mock_puller._asset.data_folder / "abc123digest" / "fasta" / asset_digest
        asset_dir.mkdir(parents=True)
        existing_row = MagicMock()
        existing_row.path = str(asset_dir)
        mock_puller._asset.get_by_digest.return_value = existing_row
        return mock_puller

    def test_the_callers_confirmer_decides(self, puller_at_the_prompt):
        asked = []

        def refuse(message):
            asked.append(message)
            return False

        with patch("refgenie.utils.prompt.ask", side_effect=AssertionError("read stdin")):
            result = puller_at_the_prompt.pull(
                asset_group_name="fasta",
                genome=GenomeAlias("rCRSd"),
                options=PullOptions(confirm=refuse),
            )

        assert result is None, "a declined overwrite is a skip, not a pull"
        assert asked and "Replace existing" in asked[0]
