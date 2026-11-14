"""``rgc.servers``: subscriptions, clients, ``locate``, the remote catalog, remote
seek, pull size estimates, and ``GenomeManager.init_from_remote``.
"""

from unittest.mock import MagicMock, patch

import pytest

from refgenie.managers.asset import AssetManager
from refgenie.managers.sources import estimate_pull_size
from refgenie.managers.sources.genomes import collection_source
from refgenie.models import GenomeAlias, GenomeDigest
from tests.support.factories import GENOME, OMIT, fake_digest, seed_alias_row
from tests.support.remotes import mock_server_client, mocked_puller
from tests.support.server import make_server_client_world, serve_refgenie


class TestSubscriptions:
    """Subscribe/unsubscribe round-trip through the Configuration row (unit tier)."""

    def test_round_trip(self, refgenie_minimal):
        servers = refgenie_minimal.servers
        assert servers.subscriptions() == []
        servers.subscribe(["http://a.example.com", "http://b.example.com"])
        assert sorted(servers.subscriptions()) == ["http://a.example.com", "http://b.example.com"]
        servers.unsubscribe(["http://a.example.com"])
        assert servers.subscriptions() == ["http://b.example.com"]
        servers.subscribe("http://c.example.com", reset=True)
        assert servers.subscriptions() == ["http://c.example.com"]


class TestSubscriptionNormalization:
    """Server URLs compare after normalization; lookups never build a client (unit tier)."""

    def test_subscribe_stores_normalized(self, refgenie_minimal):
        servers = refgenie_minimal.servers
        servers.subscribe(["http://a.example/"])
        assert servers.subscriptions() == ["http://a.example"]

    def test_find_subscription(self, refgenie_minimal):
        servers = refgenie_minimal.servers
        servers.subscribe(["http://a.example"])
        assert servers.find_subscription("http://a.example") == "http://a.example"
        assert servers.find_subscription("http://a.example/") == "http://a.example"
        assert servers.find_subscription("  http://a.example//  ") == "http://a.example"
        assert servers.find_subscription("http://b.example") is None
        assert servers.clients == {}

    def test_unsubscribe_with_trailing_slash(self, refgenie_minimal):
        servers = refgenie_minimal.servers
        servers.subscribe(["http://a.example", "http://b.example"])
        servers.unsubscribe(["http://a.example/"])
        assert servers.subscriptions() == ["http://b.example"]

    def test_unreachable_server_is_skipped_not_raised(self, refgenie_minimal):
        """Client construction fetches /openapi.json; a down server is one skipped server."""
        servers = refgenie_minimal.servers
        servers.subscribe(["http://down.example"])
        with patch(
            "refgenie.managers.sources.servers.RefgenieserverClient",
            side_effect=ConnectionError("refused"),
        ):
            assert servers.list_genomes() == []
            assert servers.list_assets_for_genome(fake_digest("g")) == []


class TestClients:
    """One client per URL, created on first use (unit tier)."""

    def test_client_is_cached(self, refgenie_minimal):
        servers = refgenie_minimal.servers
        # The real client fetches the server's OpenAPI spec when built.
        with patch(
            "refgenie.managers.sources.servers.RefgenieserverClient",
            side_effect=lambda url: MagicMock(url=url),
        ) as make_client:
            first = servers.client("http://a.example.com")
            assert servers.client("http://a.example.com") is first
        make_client.assert_called_once_with("http://a.example.com")
        assert servers.clients == {"http://a.example.com": first}


class TestLocateUnit:
    """``rgc.servers.locate`` over mocked clients (unit tier)."""

    def test_the_first_server_that_answers_wins(self, refgenie_minimal):
        down, good = MagicMock(), MagicMock()
        down.alias.resolve.side_effect = OSError("down")
        good.alias.resolve.return_value = fake_digest("a")
        good.genome.get.return_value = {"description": "hg"}
        servers = refgenie_minimal.servers
        clients = {"http://down": down, "http://good": good}
        with patch.object(servers, "client", side_effect=clients.__getitem__):
            found = servers.locate(GenomeAlias("hg38"), list(clients))
        assert found == ("http://good", fake_digest("a"), {"description": "hg"})

    def test_a_server_that_cannot_be_reached_is_skipped(self, refgenie_minimal):
        servers = refgenie_minimal.servers
        with patch.object(servers, "client", side_effect=ConnectionError("refused")):
            assert servers.locate(GenomeDigest(fake_digest("a")), ["http://down"]) is None

    def test_alias_names_is_empty_when_the_listing_fails(self, refgenie_minimal):
        client = MagicMock()
        client.alias.list.side_effect = OSError("down")
        servers = refgenie_minimal.servers
        with patch.object(servers, "client", return_value=client):
            assert servers.alias_names(GenomeDigest(fake_digest("a")), "http://a") == []


class TestCollectionSource:
    """``collection_source``: where a located genome's collection is read (unit tier)."""

    def test_reads_the_store_named_by_the_genome_record(self):
        """A server serves several stores but advertises only one; a genome in
        another store is read from the store its record names."""
        store_url = "https://store.test/vgp/"
        source = MagicMock()
        with patch("refgenie.managers.sources.genomes.make_source", return_value=source) as make:
            assert collection_source("http://server", {"refget_store_url": store_url}) is source
        make.assert_called_once_with(store_url)

    def test_with_no_store_url_it_uses_the_server(self):
        source = MagicMock()
        with patch("refgenie.managers.sources.genomes.make_source", return_value=source) as make:
            assert collection_source("http://server", {"refget_store_url": None}) is source
        make.assert_called_once_with("http://server")

    def test_a_store_that_fails_to_open_falls_back_to_the_server(self):
        source = MagicMock()
        with patch(
            "refgenie.managers.sources.genomes.make_source",
            side_effect=[ConnectionError("no store"), source],
        ) as make:
            assert collection_source("http://server", {"refget_store_url": "https://s/"}) is source
        assert [c.args[0] for c in make.call_args_list] == ["https://s/", "http://server"]

    def test_none_when_nothing_opens(self):
        with patch(
            "refgenie.managers.sources.genomes.make_source", side_effect=ConnectionError("no")
        ):
            assert collection_source("http://server", {}) is None


class TestEstimatePullSize:
    """Size estimation for bulk pulls."""

    @pytest.mark.parametrize(
        "asset_list, expected",
        [
            ([{"archive_size": 1000}, {"archive_size": 2000}, {"archive_size": 3000}], (6000, 3)),
            ([{"archive_size": 1000}, {"archive_size": None}, {"archive_size": 2000}], (3000, 3)),
            ([], (0, 0)),
        ],
        ids=["known-sizes-sum", "unknown-size-counted-not-summed", "empty-list"],
    )
    def test_estimate_pull_size(self, asset_list, expected):
        assert estimate_pull_size(asset_list) == expected


# ===========================================================================
# list-remote: querying subscribed servers for available assets
# ===========================================================================


class TestArchiveSizeChoice:
    def test_file_record_listed_first_does_not_hide_the_archive_size(self, refgenie_minimal):
        servers = refgenie_minimal.servers
        servers.subscribe(["http://a.example"])
        client = MagicMock()
        client.asset_group.list.return_value = [{"id": 1, "name": "fasta"}]
        client.asset.list.return_value = [{"digest": "ad", "name": "default"}]
        client.staged_asset.list.return_value = [
            {"mode": "file", "tarball_size": None, "tarball_digest": None},
            {"mode": "archive", "tarball_size": 868, "tarball_digest": "td"},
        ]
        with patch.object(servers, "client", return_value=client):
            [row] = servers.list_assets_for_genome(fake_digest("g"))
        assert row["archive_size"] == 868
        assert row["archive_digest"] == "td"


class TestListAssets:
    """Test list-remote (querying remote servers for available assets).

    Component tier: these build real genome folders, asset files and archives
    on disk alongside SQLite. Deselected from the bare `pytest` inner loop;
    run with `pytest -m component`.
    """

    pytestmark = pytest.mark.component

    @pytest.fixture
    def servers_world(self, engine, tmp_path, fixtures_path):
        """A client refgenie subscribed to the REAL server app for list-remote.

        Wraps the server world in create_app and injects a RefgenieserverClient
        bound to a live TestClient. The TestClient context is held open for the
        duration of the test, so remote calls made by the test hit the real
        server routers. Yields (client_rg, mock_server_url).

        Unlike the pull path, the server asset is never staged and the client
        never registers the fasta definitions -- list-remote reads catalog
        metadata only. The ``subscribe`` call is the genuine difference.
        """
        server_rg, client_rg = make_server_client_world(
            engine, tmp_path, fixtures_path, stage=False, register_client_fasta=False
        )

        # A second name on identical content: one asset row, two names (#295).
        server_rg.build.run(
            recipe_name="fasta", genome_alias="rCRSd", asset_group_name="fasta", asset_name="v2"
        )

        mock_server_url = "http://mock-refgenie-server"
        client_rg.servers.subscribe(mock_server_url)

        with serve_refgenie(client_rg, server_rg, mock_server_url):
            yield client_rg, mock_server_url

    def test_list_assets_returns_assets_and_aliases(self, servers_world):
        """list-remote returns both asset and alias data from a subscribed server."""
        client_rg, mock_server_url = servers_world
        asset_data, aliases_data = client_rg.servers.list_assets()

        # Assets: rCRSd genome's fasta assets
        assert mock_server_url in asset_data
        server_assets = asset_data[mock_server_url]
        assert len(server_assets) > 0
        for genome_digest, asset_list in server_assets.items():
            assert {"fasta:default", "fasta:v2"} <= set(asset_list)

        # Aliases: mapping that includes "rCRSd"
        assert mock_server_url in aliases_data
        server_aliases = aliases_data[mock_server_url]
        alias_values = " ".join(a for names in server_aliases.values() for a in names)
        assert "rCRSd" in alias_values

    def test_assets_table_remote(self, servers_world):
        """servers.assets_table() returns Rich tables."""
        from rich.table import Table

        client_rg, _ = servers_world
        tables = client_rg.servers.assets_table()

        assert isinstance(tables, list)
        assert len(tables) > 0
        for table in tables:
            assert isinstance(table, Table)


def test_list_assets_no_subscriptions_returns_empty(refgenie_minimal):
    """list_assets with no subscriptions returns empty dicts."""
    asset_data, aliases_data = refgenie_minimal.servers.list_assets()
    assert asset_data == {}
    assert aliases_data == {}


# --- seekr (servers.seek) file-mode URL construction ---------------------------


def _seekr_client(serving_modes, **kwargs):
    """A mock RefgenieserverClient for seekr testing (no genome_digest on the group)."""
    return mock_server_client(
        asset_group_name="fasta",
        genome_digest=OMIT,
        asset_name="test",
        asset_digest="asset_digest_001",
        serving_modes=serving_modes,
        seek_keys=[
            {"name": "fasta", "value": "rCRSd.fa", "type": "file"},
            {"name": "fai", "value": "rCRSd.fa.fai", "type": "file"},
        ],
        **kwargs,
    )


class TestSeekrFileMode:
    """seekr (seek remote) with file-mode assets, against the session-scoped
    built FASTA and a mocked server client (component tier)."""

    pytestmark = pytest.mark.component

    @pytest.mark.parametrize(
        "seek_key, expect_in_result",
        [
            (None, None),  # default seek key -> just the file endpoint
            ("fai", ".fai"),  # named seek key -> URL references that seek key's file
        ],
    )
    def test_seekr_returns_file_url_for_file_mode_asset(
        self, refgenie_session, seek_key, expect_in_result
    ):
        """seekr returns a file-endpoint URL for file-mode assets, using the requested
        seek key's file path when one is given."""
        r = refgenie_session
        mock_client = _seekr_client(serving_modes=["file"])

        kwargs = {} if seek_key is None else {"seek_key": seek_key}
        with mocked_puller(
            r,
            mock_client,
            mock_genome=False,
            mock_download_modes=False,
            mock_asset_writes=False,
        ):
            result = r.servers.seek(
                genome_digest=r.alias.resolve("rCRSd"),
                asset_group_name="fasta",
                asset_name="test",
                **kwargs,
            )

        assert result.startswith("http://test.example.com/v4/assets/asset_digest_001/files/")
        if expect_in_result is not None:
            assert expect_in_result in result

    def test_seekr_resolves_alias_via_server_without_local_genome(self, refgenie_session):
        """seekr answers for a genome the client does NOT have locally: the alias
        is resolved read-only against the server, and no local genome/alias row
        is created."""
        r = refgenie_session
        remote_digest = "0" * 32
        assert not r.alias.exists("notlocal")
        assert not r.genome.exists(remote_digest)

        mock_client = _seekr_client(serving_modes=["file"], alias_digest=remote_digest)
        with mocked_puller(
            r,
            mock_client,
            mock_genome=False,
            mock_download_modes=False,
            mock_asset_writes=False,
        ):
            result = r.servers.seek_components(r.parse_asset_registry_path("notlocal/fasta:test"))

        assert result.startswith("http://test.example.com/v4/assets/asset_digest_001/files/")
        # Read-only: nothing was written locally.
        assert not r.alias.exists("notlocal")
        assert not r.genome.exists(remote_digest)

    def test_seekr_by_digest_asks_no_server_to_resolve_an_alias(self, refgenie_session):
        """seekr by digest reaches a genome with no alias anywhere: the digest
        goes straight to the server's asset-group query."""
        r = refgenie_session
        remote_digest = fake_digest("digest-only")
        mock_client = _seekr_client(serving_modes=["file"])
        with mocked_puller(
            r,
            mock_client,
            mock_genome=False,
            mock_download_modes=False,
            mock_asset_writes=False,
        ):
            result = r.servers.seek(
                genome_digest=remote_digest, asset_group_name="fasta", asset_name="test"
            )

        assert result.startswith("http://test.example.com/v4/assets/asset_digest_001/files/")
        mock_client.get.assert_not_called()
        query = mock_client.asset_group.list.call_args.kwargs["params"]
        assert query["genome_digest"] == remote_digest

    def test_seekr_resolves_default_asset_from_server(self, refgenie_session):
        """With the asset name omitted, seekr picks the server's is_default asset
        for a non-local genome."""
        r = refgenie_session
        remote_digest = "1" * 32
        mock_client = _seekr_client(
            serving_modes=["file"], alias_digest=remote_digest, is_default=True
        )
        with mocked_puller(
            r,
            mock_client,
            mock_genome=False,
            mock_download_modes=False,
            mock_asset_writes=False,
        ):
            result = r.servers.seek_components(r.parse_asset_registry_path("notlocal/fasta"))
        assert result.startswith("http://test.example.com/v4/assets/asset_digest_001/files/")

    def test_seekr_errors_for_archive_only_asset(self, refgenie_session):
        """seekr raises ValueError for archive-only assets."""
        r = refgenie_session
        mock_client = _seekr_client(serving_modes=["archive"])

        with mocked_puller(
            r,
            mock_client,
            mock_genome=False,
            mock_download_modes=False,
            mock_asset_writes=False,
        ):
            with pytest.raises(ValueError, match="does not support file-level access"):
                r.servers.seek(
                    genome_digest=r.alias.resolve("rCRSd"),
                    asset_group_name="fasta",
                    asset_name="test",
                )


@pytest.fixture
def remote_world(engine, tmp_path, fixtures_path):
    """A bare client subscribed to the real server app (component tier).

    ``make_source`` is made to fail, so no store opens and a registered genome
    is added bare, with no collection in the client's store.
    """
    server_rg, client_rg = make_server_client_world(
        engine, tmp_path, fixtures_path, stage=False, register_client_fasta=False
    )
    url = "http://mock-refgenie-server"
    client_rg.servers.subscribe(url)
    with (
        serve_refgenie(client_rg, server_rg, url),
        patch("refgenie.managers.sources.genomes.make_source", side_effect=ConnectionError("no")),
    ):
        yield client_rg, server_rg, url


class TestLocate:
    """``rgc.servers.locate`` against the real server app (component tier)."""

    pytestmark = pytest.mark.component

    def test_an_alias_or_a_digest_the_server_knows(self, remote_world):
        client_rg, server_rg, url = remote_world
        digest = server_rg.alias.resolve(GENOME)
        for genome in (GenomeAlias(GENOME), GenomeDigest(digest)):
            found_url, found_digest, record = client_rg.servers.locate(genome)
            assert (found_url, found_digest) == (url, digest)
            assert record["digest"] == digest
        # Read-only: nothing was written here.
        assert not client_rg.genome.exists(digest)

    def test_an_unknown_alias_or_digest_is_none(self, remote_world):
        client_rg, _, _ = remote_world
        assert client_rg.servers.locate(GenomeAlias("no_such_genome")) is None
        assert client_rg.servers.locate(GenomeDigest(fake_digest("nope"))) is None

    def test_a_plain_string_is_refused(self, remote_world):
        client_rg, _, _ = remote_world
        with pytest.raises(TypeError, match="GenomeAlias or a GenomeDigest"):
            client_rg.servers.locate(GENOME)


class TestInitFromRemote:
    """``GenomeManager.init_from_remote`` against the real server app (component tier)."""

    pytestmark = pytest.mark.component

    def test_registers_genome_and_alias(self, remote_world, caplog):
        """The client's store lacks rCRSd's collection, so require_collection
        undoes what the call created and suggests a fasta pull (#299)."""
        client_rg, server_rg, _ = remote_world
        digest = server_rg.alias.resolve(GENOME)
        alias = GenomeAlias(GENOME)
        assert client_rg.genome.init_from_remote(alias, require_collection=True) is None
        assert not client_rg.genome.exists(digest)
        assert "refgenie pull rCRSd/fasta" in caplog.text

        assert client_rg.genome.init_from_remote(alias) == digest
        assert client_rg.alias.resolve(alias) == digest
        assert client_rg.genome.exists(digest)

        client_rg.alias.remove(alias)  # genome exists, alias will be new
        assert client_rg.genome.init_from_remote(alias, require_collection=True) is None
        assert client_rg.genome.exists(digest)
        assert not client_rg.alias.exists(alias)

    def test_a_digest_is_registered_under_the_servers_first_usable_alias(self, remote_world):
        client_rg, server_rg, _ = remote_world
        digest = server_rg.alias.resolve(GENOME)
        seed_alias_row(server_rg, "Homo_sapiens/NCBI/GRCh38", digest)
        assert client_rg.genome.init_from_remote(GenomeDigest(digest)) == digest
        assert [str(a) for a in client_rg.alias.get_for_genome(digest)] == [GENOME]

    def test_a_digest_with_only_a_slash_alias_is_registered_under_its_digest(self, remote_world):
        client_rg, server_rg, _ = remote_world
        digest = server_rg.alias.resolve(GENOME)
        seed_alias_row(server_rg, "Homo_sapiens/NCBI/GRCh38", digest)
        server_rg.alias.remove(GENOME)
        assert client_rg.genome.init_from_remote(GenomeDigest(digest)) == digest
        assert [str(a) for a in client_rg.alias.get_for_genome(digest)] == [digest]

    def test_a_digest_already_here_asks_no_server(self, remote_world):
        client_rg, _, _ = remote_world
        digest = GenomeDigest("d" * 32)
        client_rg.genome.add(digest, "local", [])
        with patch.object(client_rg.servers, "locate", side_effect=AssertionError("asked")):
            assert client_rg.genome.init_from_remote(digest) == digest

    def test_existing_alias_is_a_no_op(self, remote_world):
        client_rg, _, _ = remote_world
        client_rg.genome.add("d" * 32, "local", [GENOME])
        with patch.object(client_rg.genome, "ensure_from_remote") as ensure:
            assert client_rg.genome.init_from_remote(GenomeAlias(GENOME)) == "d" * 32
        ensure.assert_not_called()

    def test_unknown_alias_returns_none(self, remote_world):
        client_rg, _, _ = remote_world
        assert client_rg.genome.init_from_remote(GenomeAlias("no_such_genome")) is None
        assert not client_rg.alias.exists("no_such_genome")


def test_asset_manager_needs_no_source_or_recipe_manager(refgenie_minimal):
    """AssetManager builds from records-only dependencies; it creates no builder
    or puller, so it takes no source or recipe manager."""
    r = refgenie_minimal
    asset = AssetManager(
        database_engine=r.database_engine,
        genome_folder=r.genome_folder,
        alias_folder=r.alias_folder,
        alias_manager=r.alias,
        genome_manager=r.genome,
        asset_class_manager=r.asset_class,
    )
    assert list(asset.list_assets()) == []
    assert asset.links is not None
