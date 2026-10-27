"""
Tests for ``rgc.servers`` (``refgenie.managers.sources.servers.ServerManager``):
subscriptions, clients, the remote catalog, remote seek, and pull-size
estimates. Also the genome bootstrap that uses it
(``GenomeManager.init_from_remote``), and the proof that ``AssetManager`` no
longer needs source or recipe managers.

Remote servers are the REAL server app served in-process through
``serve_refgenie``; seekr file-mode tests use a mocked client.
"""

from unittest.mock import MagicMock, patch

import pytest

from refgenie.managers.asset import AssetManager
from refgenie.managers.sources import ServerManager, estimate_pull_size
from tests.helpers import (
    OMIT,
    fake_digest,
    make_server_client_world,
    mock_server_client,
    mocked_puller,
    serve_refgenie,
)


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

    def test_injected_clients_are_validated(self, refgenie_minimal):
        with pytest.raises(ValueError, match="ServerClient Protocol"):
            ServerManager(
                refgenie_minimal.database_engine,
                alias_manager=refgenie_minimal.alias,
                seek_keys=refgenie_minimal.asset.seek_key,
                clients={"http://a.example.com": object()},
            )


class TestGenomeDescription:
    """The first server that describes a genome wins (unit tier)."""

    def test_skips_failing_and_empty_servers(self, refgenie_minimal):
        clients = {
            "http://down.example.com": MagicMock(get=MagicMock(side_effect=OSError("down"))),
            "http://empty.example.com": MagicMock(get=MagicMock(return_value={"description": ""})),
            "http://good.example.com": MagicMock(get=MagicMock(return_value={"description": "hg"})),
        }
        servers = refgenie_minimal.servers
        with patch.object(servers, "client", side_effect=clients.__getitem__):
            assert servers.genome_description("abc", list(clients)) == "hg"

    def test_empty_when_no_server_describes_it(self, refgenie_minimal):
        servers = refgenie_minimal.servers
        client = MagicMock(get=MagicMock(return_value={}))
        with patch.object(servers, "client", return_value=client):
            assert servers.genome_description("abc", ["http://a.example.com"]) == ""


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
            assert len(asset_list) > 0
            assert any("fasta:" in a for a in asset_list)

        # Aliases: mapping that includes "rCRSd"
        assert mock_server_url in aliases_data
        server_aliases = aliases_data[mock_server_url]
        alias_values = " ".join(server_aliases.values())
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
    """A mock ServerClient for seekr testing (no genome_digest on the group)."""
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
        query = mock_client.get_asset_groups.call_args.kwargs["params"]
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


class TestFindCollection:
    """``rgc.servers.find_collection``, behind ``refgenie id --remote`` (unit tier)."""

    def test_returns_none_with_no_servers(self, refgenie_minimal):
        assert refgenie_minimal.servers.find_collection(fake_digest("a"), server_urls=[]) is None

    def test_returns_none_on_exception(self, refgenie_minimal):
        with patch(
            "refgenie.managers.sources.servers.make_source",
            side_effect=ConnectionError("unreachable"),
        ) as mock_make_source:
            result = refgenie_minimal.servers.find_collection(
                fake_digest("a"), server_urls=["http://fake.server"]
            )
        assert result is None
        mock_make_source.assert_called_once_with("http://fake.server")

    def test_returns_none_when_source_lacks_collection(self, refgenie_minimal):
        source = MagicMock()
        source.verify_collection.return_value = None
        with patch("refgenie.managers.sources.servers.make_source", return_value=source):
            result = refgenie_minimal.servers.find_collection(
                fake_digest("a"), server_urls=["http://fake.server"]
            )
        assert result is None
        source.verify_collection.assert_called_once_with(fake_digest("a"))

    def test_returns_the_first_servers_collection(self, refgenie_minimal):
        missing, found = MagicMock(), MagicMock()
        missing.verify_collection.return_value = None
        found.verify_collection.return_value = {"names": ["chrM"]}
        with patch(
            "refgenie.managers.sources.servers.make_source",
            side_effect=lambda url: {"http://a": missing, "http://b": found}[url],
        ):
            result = refgenie_minimal.servers.find_collection(
                fake_digest("a"), server_urls=["http://a", "http://b"]
            )
        assert result == {"names": ["chrM"]}


class TestResolveAlias:
    """``rgc.servers.resolve_alias``: the genome source first, then the v4 client (unit tier)."""

    def test_source_answers_first(self, refgenie_minimal):
        source = MagicMock()
        source.resolve_alias.return_value = fake_digest("a")
        client = MagicMock()
        with (
            patch("refgenie.managers.sources.servers.make_source", return_value=source),
            patch.object(refgenie_minimal.servers, "client", return_value=client),
        ):
            assert refgenie_minimal.servers.resolve_alias("hg38", ["http://a"]) == fake_digest("a")
        client.get.assert_not_called()

    def test_client_is_the_fallback(self, refgenie_minimal):
        client = MagicMock()
        client.get.return_value = {"digest": fake_digest("b")}
        with (
            patch(
                "refgenie.managers.sources.servers.make_source",
                side_effect=ConnectionError("no"),
            ),
            patch.object(refgenie_minimal.servers, "client", return_value=client),
        ):
            assert refgenie_minimal.servers.resolve_alias("hg38", ["http://a"]) == fake_digest("b")

    def test_unknown_everywhere_is_none(self, refgenie_minimal):
        source = MagicMock()
        source.resolve_alias.return_value = None
        client = MagicMock()
        client.get.side_effect = ValueError("not found")
        with (
            patch("refgenie.managers.sources.servers.make_source", return_value=source),
            patch.object(refgenie_minimal.servers, "client", return_value=client),
        ):
            assert refgenie_minimal.servers.resolve_alias("hg38", ["http://a"]) is None


class TestInitFromRemote:
    """``GenomeManager.init_from_remote`` against the real server app (component tier).

    ``make_source`` is made to fail so the genome is resolved through the
    server client's alias endpoint, which the in-process app answers.
    """

    pytestmark = pytest.mark.component

    @pytest.fixture
    def world(self, engine, tmp_path, fixtures_path):
        server_rg, client_rg = make_server_client_world(
            engine, tmp_path, fixtures_path, stage=False, register_client_fasta=False
        )
        url = "http://mock-refgenie-server"
        client_rg.servers.subscribe(url)
        with (
            serve_refgenie(client_rg, server_rg, url),
            patch(
                "refgenie.managers.sources.servers.make_source", side_effect=ConnectionError("no")
            ),
        ):
            yield client_rg, server_rg

    def test_registers_genome_and_alias(self, world):
        client_rg, server_rg = world
        assert not client_rg.alias.exists("rCRSd")
        assert client_rg.genome.init_from_remote("rCRSd") is True
        assert client_rg.alias.resolve("rCRSd") == server_rg.alias.resolve("rCRSd")
        assert client_rg.genome.exists(server_rg.alias.resolve("rCRSd"))

    def test_existing_alias_is_a_no_op(self, world):
        client_rg, _ = world
        client_rg.genome.add("d" * 32, "local", ["rCRSd"])
        with patch.object(client_rg.genome, "ensure_from_remote") as ensure:
            assert client_rg.genome.init_from_remote("rCRSd") is True
        ensure.assert_not_called()

    def test_unknown_alias_returns_false(self, world):
        client_rg, _ = world
        assert client_rg.genome.init_from_remote("no_such_genome") is False
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
