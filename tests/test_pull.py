"""Tests for the refgenie pull workflow -- pulling assets from remote servers.

Organized by behavior:

* TestPullFromServer          -- full Refgenie.pull() against the REAL server app
                                 (create_app) served in-process via serve_refgenie
* TestPullFileMode            -- the _pull_file_mode file-by-file download helper
* TestPullDownloadModeSelect  -- archive/file/none serving-mode dispatch in pull()
* TestPullFailures            -- pull() raises (never returns a silent None) on
                                 unreachable/empty/ambiguous server responses
* TestArchiveVerification     -- archive checksum verification guard
* TestPullRollback            -- atomicity: no orphaned genome/alias/dirs on failure
* TestBuildPullParents        -- build_asset(pull_parents=True) delegates to pull()
* TestListRemote              -- `listr` against the real server app
* TestEstimatePullSize        -- size estimation for bulk pulls
* TestPullMultiple            -- pull_multiple failure handling

Bulk-pull argv parsing and the bulk-pull CLI exit codes live in test_cli.py; the
confirmation-prompt utilities live in test_utils.py.

Tier markers are applied PER CLASS: the server/filesystem-heavy classes are
``component``; the mock-only classes default to ``unit`` (assigned by conftest).
"""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from tests.helpers import (
    MOCK_MIRROR_URL,
    MOCK_SERVER_URL,
    OMIT,
    MockRemoteSource,
    make_engine,
    make_server_client_world,
    mock_server_client,
    mocked_puller,
    requires_server,
    serve_refgenie,
)


class TestPullFromServer:
    """Full Refgenie.pull() exercised end-to-end against the REAL server app.

    ``server_client_world`` (tests/conftest.py) serves the server through the
    REAL app (create_app) over an in-process ASGI transport -- no hand-rolled
    routes.
    """

    pytestmark = pytest.mark.component
    # The real app (server extras) is required only by this class; guard it here
    # so the mock-only classes still run without it.
    requires_server()

    def test_pull_asset_from_server(self, server_client_world):
        """Pulling a FASTA asset from a real remote server app registers it locally."""
        client_rg, _, url = server_client_world
        result = client_rg.pull(
            alias_name="rCRSd",
            asset_group_name="fasta",
            force=True,
            force_large=True,
            force_server_urls=[url],
        )

        assert result is not None, "pull() returned None - asset was not pulled"
        assert client_rg.asset.exists(
            genome_name="rCRSd", asset_group_name="fasta", asset_name="default"
        )
        assert client_rg.alias.exists("rCRSd")
        assert client_rg.genome.exists(client_rg.alias.resolve("rCRSd"))

    def test_pull_carries_the_servers_build_provenance(self, server_client_world):
        """A pulled name records the build the server recorded for it.

        Without this, ``build_digest`` is only addressable on the machine that
        built the asset, which defeats making builds addressable at all.
        """
        from tests.helpers import asset_name_rows

        client_rg, server_rg, url = server_client_world
        served = {
            row.name: row for row in asset_name_rows(server_rg) if row.build_digest is not None
        }
        assert served, "the server fixture must have built something to serve"

        client_rg.pull(
            alias_name="rCRSd",
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
            # client's, so the recipe link is deliberately not carried over.
            assert pulled[name].recipe_id is None

    def test_pull_stops_after_first_server_succeeds(self, server_client_world):
        """Pull with several subscribed servers must stop at the first that works.

        Regression: the server loop had no ``break``, so after a successful pull
        from server A the loop re-entered for server B with the asset now
        registered locally, raising AssetExistsError for a pull that had just
        succeeded.
        """
        client_rg, _, _ = server_client_world
        # force is left at its default (None) -- that is the code path where the
        # missing break turned into a user-visible AssetExistsError.
        result = client_rg.pull(
            alias_name="rCRSd",
            asset_group_name="fasta",
            force_large=True,
            force_server_urls=[MOCK_SERVER_URL, MOCK_MIRROR_URL],
        )

        assert result is not None
        assert client_rg.asset.exists(
            genome_name="rCRSd", asset_group_name="fasta", asset_name="default"
        )

    def test_pull_multiple_counts_multi_server_pull_as_success(self, server_client_world):
        """pull_multiple must book a success when several servers are subscribed.

        Regression for the missing ``break``: without it the second iteration
        raised AssetExistsError and pull_multiple recorded a failure even though
        the asset had been downloaded and registered.
        """
        client_rg, server_rg, _ = server_client_world
        urls = [MOCK_SERVER_URL, MOCK_MIRROR_URL]
        genome_digest = server_rg.alias.resolve("rCRSd")

        puller = client_rg.asset._asset_puller
        with patch.object(client_rg.sources, "get_subscriptions", return_value=urls):
            successful = puller.pull_multiple(
                asset_list=[
                    {
                        "genome_digest": genome_digest,
                        "asset_group_name": "fasta",
                        "asset_name": "default",
                    }
                ],
                force_large=True,
            )

        assert len(successful) == 1, "multi-server pull should be counted as a success"

    def test_existing_asset_raises_without_force(self, server_client_world):
        """Pulling an already-present asset raises without force, replaces with it."""
        from refgenie.exceptions import AssetExistsError

        client_rg, _, url = server_client_world
        kwargs = dict(alias_name="rCRSd", asset_group_name="fasta", force_large=True,
                      force_server_urls=[url])

        assert client_rg.pull(force=True, **kwargs) is not None

        with pytest.raises(AssetExistsError):
            client_rg.pull(force=False, **kwargs)

        # force=True over the existing asset re-downloads and replaces it.
        assert client_rg.pull(force=True, **kwargs) is not None
        assert client_rg.asset.exists(
            genome_name="rCRSd", asset_group_name="fasta", asset_name="default"
        )


# ===========================================================================
# File-mode download + serving-mode dispatch
# ===========================================================================


_MOCK_ALIASES = [{"name": "test_genome", "genome_digest": "genome_digest_001"}]


class TestPullFileMode:
    """Tests for the _pull_file_mode file-by-file download helper."""

    pytestmark = pytest.mark.component

    def test_downloads_files_into_nested_dirs(self, refgenie_fs, tmp_path):
        """_pull_file_mode downloads each file from the list, creating subdirs."""
        from refgenie.managers.asset.puller import PullTransaction

        puller = refgenie_fs.asset._asset_puller
        client = MagicMock()
        client.get_asset_file_list.return_value = ["genome.gtf", "subdir/annotation.bed"]

        def mock_download(asset_digest, file_path, output_path):
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text(f"content of {file_path}")
            return output_path

        client.download_file.side_effect = mock_download

        asset_dir = tmp_path / "test_asset"
        puller._pull_file_mode(
            client=client,
            asset_digest="abc123",
            asset_dir=asset_dir,
            bundle_name="test/fasta:default",
            txn=PullTransaction(puller),
        )

        assert client.download_file.call_count == 2
        assert (asset_dir / "genome.gtf").exists()
        assert (asset_dir / "subdir" / "annotation.bed").exists()

    def test_empty_file_list_raises(self, refgenie_fs, tmp_path):
        """_pull_file_mode raises RuntimeError (never downloads) on an empty list."""
        from refgenie.managers.asset.puller import PullTransaction

        puller = refgenie_fs.asset._asset_puller
        client = MagicMock()
        client.get_asset_file_list.return_value = []

        with pytest.raises(RuntimeError, match="No files found"):
            puller._pull_file_mode(
                client=client,
                asset_digest="abc123",
                asset_dir=tmp_path / "test_asset",
                bundle_name="test/fasta:default",
                txn=PullTransaction(puller),
            )
        client.download_file.assert_not_called()

    def test_rollback_removes_created_dir_on_failure(self, refgenie_fs, tmp_path):
        """On a mid-download failure, PullTransaction removes the created directory."""
        from refgenie.managers.asset.puller import PullTransaction

        puller = refgenie_fs.asset._asset_puller
        client = MagicMock()
        client.get_asset_file_list.return_value = ["file1.txt", "file2.txt"]

        call_count = [0]

        def mock_download(asset_digest, file_path, output_path):
            call_count[0] += 1
            if call_count[0] == 1:
                output_path.parent.mkdir(parents=True, exist_ok=True)
                output_path.write_text("content")
                return output_path
            raise ConnectionError("Download failed")

        client.download_file.side_effect = mock_download

        asset_dir = tmp_path / "rollback_test"
        txn = PullTransaction(puller)
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


# Archive-mode staged record (byte-digest + size present).
_ARCHIVE_REC = {
    "id": 1,
    "asset_digest": "abc123digest",
    "mode": "archive",
    "tarball_size": 1000,
    "tarball_digest": "tardigest",
}
# File-mode staged record.
_FILE_REC = {
    "id": 1,
    "asset_digest": "abc123digest",
    "mode": "file",
    "directory_contents": ["genome.fa"],
}


class TestPullDownloadModeSelect:
    """pull() must dispatch to the right download helper for the serving modes."""

    pytestmark = pytest.mark.component

    @pytest.mark.parametrize(
        "serving_modes, staged, file_list, expected",
        [
            # Both modes offered, archive is staged -> prefer archive.
            (["file", "archive"], [_ARCHIVE_REC], [], "archive"),
            # Only the file mode offered and staged -> file.
            (["file"], [_FILE_REC], ["genome.fa"], "file"),
            # serving_modes omitted entirely -> fall back to what is staged (file).
            (None, [_FILE_REC], ["genome.fa"], "file"),
            # Both declared but only a file record staged -> file.
            (["file", "archive"], [_FILE_REC], ["genome.fa"], "file"),
        ],
        ids=["both_prefer_archive", "only_file", "modes_absent_fallback", "archive_declared_not_staged"],
    )
    def test_selects_download_mode(
        self, refgenie_fs, serving_modes, staged, file_list, expected
    ):
        """Two server URLs are passed so assert_called_once also guards the loop
        ``break`` -- without it the second server would be downloaded from too.
        """
        mock_client = mock_server_client(
            serving_modes=OMIT if serving_modes is None else serving_modes,
            staged_items=staged,
            file_list=file_list,
            aliases=_MOCK_ALIASES,
        )

        with mocked_puller(refgenie_fs, mock_client) as m:
            m.puller.pull(
                alias_name="rCRSd",
                asset_group_name="fasta",
                force_server_urls=["http://test.example.com", "http://mirror.example.com"],
            )

        if expected == "archive":
            m.archive.assert_called_once()
            m.file.assert_not_called()
        else:
            m.file.assert_called_once()
            m.archive.assert_not_called()

    def test_none_mode_raises_pull_failed(self, refgenie_fs):
        """serving_modes=['none'] (metadata-only) raises PullFailedError."""
        from refgenie.exceptions import PullFailedError

        mock_client = mock_server_client(serving_modes=["none"], aliases=_MOCK_ALIASES)

        with mocked_puller(
            refgenie_fs, mock_client, mock_download_modes=False, mock_asset_writes=False
        ) as m:
            with pytest.raises(PullFailedError, match="metadata-only"):
                m.puller.pull(
                    alias_name="rCRSd",
                    asset_group_name="fasta",
                    force_server_urls=["http://test.example.com"],
                )


# ===========================================================================
# Failure surfacing (mock-only, unit tier)
# ===========================================================================


@pytest.fixture
def mock_puller(tmp_path):
    """An AssetPuller with fully mocked collaborators and one server subscription.

    Distinct from ``tests.helpers.mocked_puller``: this *constructs* a real
    AssetPuller over MagicMock collaborators, rather than patching an existing
    puller's methods. Both strategies are legitimate; do not merge them.
    """
    from refgenie.exceptions import MissingAssetGroupError
    from refgenie.managers.asset.puller import AssetPuller

    engine = make_engine()
    source_manager = MagicMock()
    asset_manager = MagicMock()
    alias_manager = MagicMock()
    genome_manager = MagicMock()
    asset_relations = MagicMock()

    puller = AssetPuller(
        database_engine=engine,
        genome_folder=tmp_path / "genomes",
        alias_folder=tmp_path / "aliases",
        source_manager=source_manager,
        asset_manager=asset_manager,
        alias_manager=alias_manager,
        genome_manager=genome_manager,
        asset_relations=asset_relations,
    )
    source_manager.get_subscriptions.return_value = ["http://test-server:5000"]
    alias_manager.resolve.return_value = "abc123digest"
    alias_manager.exists.return_value = True
    asset_manager.exists.return_value = False
    asset_manager.get_default.side_effect = MissingAssetGroupError("genome", "fasta")
    return puller


class TestPullFailures:
    """pull() must raise or log-and-continue on bad server responses, never
    silently return None where a real error occurred."""

    @pytest.mark.parametrize(
        "exc",
        [ValueError("Missing operation ID"), ConnectionError("Connection refused")],
        ids=["value_error", "connection_error"],
    )
    def test_unreachable_asset_groups_logs_and_returns_none(self, mock_puller, caplog, exc):
        """When get_asset_groups raises, the single server is exhausted; pull logs
        the failure and returns None (no other server to fall through to)."""
        client = MagicMock()
        client.get_asset_groups.side_effect = exc
        mock_puller._sources.get_server_client.return_value = client

        result = mock_puller.pull(alias_name="testgenome", asset_group_name="fasta")
        assert result is None
        assert "Failed to query asset groups from http://test-server:5000" in caplog.text

    def test_zero_asset_groups_raises(self, mock_puller):
        """Zero asset groups returned raises PullFailedError."""
        from refgenie.exceptions import PullFailedError

        client = MagicMock()
        client.get_asset_groups.return_value = []
        mock_puller._sources.get_server_client.return_value = client

        with pytest.raises(PullFailedError, match="Expected one asset group, got 0"):
            mock_puller.pull(alias_name="testgenome", asset_group_name="fasta")

    def test_multiple_assets_no_name_raises_helpful(self, mock_puller):
        """Multiple assets with no name specified raises with a 'specify one' hint."""
        from refgenie.exceptions import PullFailedError

        client = MagicMock()
        client.get_asset_groups.return_value = [
            {"id": "group1", "name": "fasta", "genome_digest": "abc123digest"}
        ]
        client.get_assets.return_value = [
            {"id": "asset1", "name": "default", "digest": "d1"},
            {"id": "asset2", "name": "v2", "digest": "d2"},
        ]
        mock_puller._sources.get_server_client.return_value = client

        with pytest.raises(PullFailedError, match="Multiple assets found.*Specify one with"):
            mock_puller.pull(alias_name="testgenome", asset_group_name="fasta", asset_name=None)

    def test_genome_resolution_failure_logs_warning(self, mock_puller, caplog):
        """When the genome alias cannot be resolved, log a warning with the URL."""
        from refgenie.exceptions import MissingAliasError
        from refgenie.managers.asset.puller import GenomeCreationResult

        mock_puller._alias_manager.resolve.side_effect = MissingAliasError("testgenome")
        with patch.object(
            mock_puller, "_ensure_genome_exists", return_value=GenomeCreationResult(success=False)
        ):
            result = mock_puller.pull(alias_name="testgenome", asset_group_name="fasta")

        assert result is None
        assert "Could not resolve genome 'testgenome' via server" in caplog.text
        assert "http://test-server:5000" in caplog.text

    def test_no_servers_raises(self, mock_puller):
        """No server subscriptions raises PullFailedError."""
        from refgenie.exceptions import PullFailedError

        mock_puller._sources.get_subscriptions.return_value = []
        with pytest.raises(PullFailedError, match="No server subscriptions found"):
            mock_puller.pull(alias_name="testgenome", asset_group_name="fasta")


class TestArchiveVerification:
    """Archive verification must use the tarball byte-digest, or say it did not run."""

    def test_warns_when_tarball_digest_absent(self, mock_puller, tmp_path, caplog):
        """A staged record without ``tarball_digest`` must warn, not silently skip.

        It must NOT fall back to the record's ``digest`` field -- that is the asset
        *identity* digest, and comparing it against the tarball's byte digest would
        produce a bogus 'checksum mismatch'.
        """
        import logging

        from refgenie.managers.asset.puller import PullTransaction

        client = MagicMock()

        def fake_download(operation_id, output_path, url_format_params, name):
            Path(output_path).parent.mkdir(parents=True, exist_ok=True)
            Path(output_path).write_bytes(b"tarball bytes")

        client.download_with_progress.side_effect = fake_download

        staged = {"mode": "archive", "digest": "asset-identity-digest", "size": 999}
        with caplog.at_level(logging.WARNING):
            with patch("refgenie.managers.asset.puller.untar"):
                mock_puller._pull_archive_mode(
                    client=client,
                    staged_asset_metadata=staged,
                    asset_digest="asset-identity-digest",
                    genome_digest="g1",
                    asset_group_name="fasta",
                    asset_name="default",
                    asset_dir=tmp_path / "data" / "g1" / "fasta" / "asset-identity-digest",
                    bundle_name="g1/fasta:default",
                    txn=PullTransaction(mock_puller),
                    force_large=True,
                    size_cutoff=None,
                    sigint_handler=None,
                )

        assert "could not be verified" in caplog.text.lower()


class TestPullSigintHandling:
    """The SIGINT handler installed around the archive download must be
    restored afterwards, not left permanently rebound to a stale tarball path.
    """

    def _run(self, mock_puller, tmp_path, *, download_side_effect=None):
        from refgenie.managers.asset.puller import PullTransaction

        client = MagicMock()
        if download_side_effect is not None:
            client.download_with_progress.side_effect = download_side_effect
        else:

            def fake_download(operation_id, output_path, url_format_params, name):
                Path(output_path).parent.mkdir(parents=True, exist_ok=True)
                Path(output_path).write_bytes(b"tarball bytes")

            client.download_with_progress.side_effect = fake_download

        staged = {"mode": "archive", "tarball_digest": None, "tarball_size": 999}
        with patch("refgenie.managers.asset.puller.untar"):
            mock_puller._pull_archive_mode(
                client=client,
                staged_asset_metadata=staged,
                asset_digest="asset-digest",
                genome_digest="g1",
                asset_group_name="fasta",
                asset_name="default",
                asset_dir=tmp_path / "data" / "g1" / "fasta" / "asset-digest",
                bundle_name="g1/fasta:default",
                txn=PullTransaction(mock_puller),
                force_large=True,
                size_cutoff=None,
                sigint_handler=None,
            )

    def test_sigint_handler_restored_after_successful_download(self, mock_puller, tmp_path):
        """Regression: the handler installed for the download was never restored,
        so the process's SIGINT stayed bound to the tarball path forever after."""
        import signal

        original_handler = signal.getsignal(signal.SIGINT)
        try:
            self._run(mock_puller, tmp_path)
            assert signal.getsignal(signal.SIGINT) is original_handler
        finally:
            signal.signal(signal.SIGINT, original_handler)

    def test_sigint_handler_restored_after_failed_download(self, mock_puller, tmp_path):
        """The handler must be restored even when the download raises -- the
        try/finally must cover the exception path, not just the happy path."""
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
# Atomicity / rollback (mock-only, unit tier)
# ===========================================================================


# The one collection every rollback test's remote source advertises.
_ROLLBACK_COLLECTION = {"names": ["chr1"], "lengths": [1000], "sequences": ["seq1"]}
_ROLLBACK_URL = "http://example.com"


def _rollback_client(asset_group_name, genome_digest, staged, download_side_effect=None):
    """A mock server client for the rollback tests.

    ``serving_modes=OMIT`` is load-bearing: the server declares none, so the
    puller falls back to what is actually staged.
    """
    return mock_server_client(
        server_url=_ROLLBACK_URL,
        asset_group_name=asset_group_name,
        genome_digest=genome_digest,
        asset_digest="def456",
        serving_modes=OMIT,
        staged_items=staged,
        relations={"parents": []},
        download_side_effect=download_side_effect,
    )


def _rollback_source(genome_digest, alias_name):
    """A real-protocol remote source that resolves ``alias_name`` locally."""
    return MockRemoteSource(
        collections={genome_digest: _ROLLBACK_COLLECTION},
        store_url=f"{_ROLLBACK_URL}/store",
        aliases={alias_name: genome_digest},
    )


def _rollback_patches(rg, client, genome_digest, alias_name, *, remote=True):
    """The shared rollback harness: real genome writes, mocked server + source."""
    return mocked_puller(
        rg,
        client,
        subscriptions=[_ROLLBACK_URL],
        remote_source=_rollback_source(genome_digest, alias_name) if remote else None,
        mock_genome=False,
        mock_download_modes=False,
        mock_asset_writes=False,
    )


class TestPullRollback:
    """A failed pull must leave no orphaned genome/alias/directory behind, and
    must not touch state that existed before the pull."""

    def test_rollback_on_download_failure(self, refgenie_minimal):
        """A download failure rolls back the newly created genome and alias."""
        rg = refgenie_minimal
        alias_name, asset_group_name, genome_digest = "test_genome", "bowtie2_index", "abc123"
        assert not rg.alias.exists(alias_name)
        assert not rg.genome.exists(genome_digest)

        client = _rollback_client(
            asset_group_name,
            genome_digest,
            staged=[{"tarball_digest": "archive123", "tarball_size": 1000, "mode": "archive"}],
            download_side_effect=Exception("Network error"),
        )

        with _rollback_patches(rg, client, genome_digest, alias_name):
            with pytest.raises(Exception, match="Network error"):
                rg.pull(alias_name=alias_name, asset_group_name=asset_group_name)

        assert not rg.alias.exists(alias_name), "Alias should be rolled back on failure"
        assert not rg.genome.exists(genome_digest), "Genome should be rolled back on failure"

    def test_rollback_on_checksum_mismatch(self, refgenie_minimal):
        """A downloaded archive whose checksum does not match rolls everything back."""
        rg = refgenie_minimal
        alias_name, asset_group_name, genome_digest = "test_genome2", "bowtie2_index", "xyz789"
        assert not rg.alias.exists(alias_name)
        assert not rg.genome.exists(genome_digest)

        def fake_download(operation_id, output_path, url_format_params, name):
            Path(output_path).parent.mkdir(parents=True, exist_ok=True)
            Path(output_path).write_bytes(b"different content that wont match checksum")

        client = _rollback_client(
            asset_group_name,
            genome_digest,
            staged=[
                {
                    "tarball_digest": "expected_checksum_that_wont_match",
                    "tarball_size": 1000,
                    "mode": "archive",
                }
            ],
            download_side_effect=fake_download,
        )

        with _rollback_patches(rg, client, genome_digest, alias_name):
            with pytest.raises(ValueError, match="checksum mismatch"):
                rg.pull(alias_name=alias_name, asset_group_name=asset_group_name)

        assert not rg.alias.exists(alias_name), "Alias should be rolled back on failure"
        assert not rg.genome.exists(genome_digest), "Genome should be rolled back on failure"

    def test_rollback_when_server_has_no_archive(self, refgenie_minimal):
        """A server with nothing staged must not leave a genome/alias/dir behind.

        Regression: ``continue``/``return None`` exit a ``with`` block normally, so
        PullTransaction.__exit__ never rolled back. The 'no staged archive' path
        therefore created a genome + alias and left them there.
        """
        from refgenie.exceptions import NoArchiveError

        rg = refgenie_minimal
        alias_name, asset_group_name, genome_digest = (
            "no_archive_genome",
            "bowtie2_index",
            "noarchive123",
        )
        assert not rg.alias.exists(alias_name)
        assert not rg.genome.exists(genome_digest)

        client = _rollback_client(asset_group_name, genome_digest, staged=[])

        with _rollback_patches(rg, client, genome_digest, alias_name):
            with pytest.raises(NoArchiveError):
                rg.pull(alias_name=alias_name, asset_group_name=asset_group_name)

        assert not rg.alias.exists(alias_name), "Alias should be rolled back"
        assert not rg.genome.exists(genome_digest), "Genome should be rolled back"
        group_dir = Path(rg.genome_folder) / "data" / genome_digest
        assert not group_dir.exists(), f"Orphaned directory left behind: {group_dir}"

    def test_rollback_preserves_genome_known_under_other_alias(self, refgenie_minimal):
        """A failed pull that added a NEW alias to a pre-existing genome must
        roll back only the alias, never the genome.

        Regression: when the alias had to be resolved remotely,
        _ensure_genome_exists reported created_genome=True even though the
        resolved digest already existed locally under another alias, so
        rollback deleted a genome the pull did not create.
        """
        rg = refgenie_minimal
        asset_group_name, genome_digest = "bowtie2_index", "shared123"
        rg.genome.add(
            digest=genome_digest, description="Pre-existing genome", alias_names=["old_alias"]
        )
        assert not rg.alias.exists("new_alias")

        client = _rollback_client(
            asset_group_name,
            genome_digest,
            staged=[{"tarball_digest": "archive123", "tarball_size": 1000, "mode": "archive"}],
            download_side_effect=Exception("Network error"),
        )

        with _rollback_patches(rg, client, genome_digest, "new_alias"):
            with pytest.raises(Exception, match="Network error"):
                rg.pull(alias_name="new_alias", asset_group_name=asset_group_name)

        assert rg.genome.exists(genome_digest), "Pre-existing genome must survive rollback"
        assert rg.alias.exists("old_alias"), "Pre-existing alias must survive rollback"
        assert not rg.alias.exists("new_alias"), "Alias created by the failed pull is rolled back"

    def test_no_rollback_when_genome_preexists(self, refgenie_minimal):
        """A pull failure must NOT delete a genome/alias that existed beforehand."""
        rg = refgenie_minimal
        alias_name, asset_group_name, genome_digest = (
            "existing_genome",
            "bowtie2_index",
            "preexisting123",
        )
        rg.genome.add(digest=genome_digest, description="Pre-existing genome", alias_names=[alias_name])
        assert rg.alias.exists(alias_name)
        assert rg.genome.exists(genome_digest)

        client = _rollback_client(
            asset_group_name,
            genome_digest,
            staged=[{"tarball_digest": "archive123", "tarball_size": 1000, "mode": "archive"}],
            download_side_effect=Exception("Network error"),
        )

        # No remote source here, deliberately: the genome and alias already exist
        # locally, so _ensure_genome_exists resolves without ever calling
        # make_source. Its four siblings patch it because they do not.
        with _rollback_patches(rg, client, genome_digest, alias_name, remote=False):
            with pytest.raises(Exception, match="Network error"):
                rg.pull(alias_name=alias_name, asset_group_name=asset_group_name)

        assert rg.alias.exists(alias_name), "Pre-existing alias should not be rolled back"
        assert rg.genome.exists(genome_digest), "Pre-existing genome should not be rolled back"


# ===========================================================================
# Bulk pull helpers (mock-only, unit)
# ===========================================================================


class TestBuildPullParents:
    """build_asset(pull_parents=True) delegates to pull() with the right arguments."""

    def test_build_asset_pull_parents_passes_genome_as_alias(self, refgenie_minimal, monkeypatch):
        """``--pull-parents`` must pull the *fasta* group for the *genome* alias.

        Regression test: build_asset once called pull() with the two arguments
        swapped, requesting the genome as an asset group.
        """
        from refgenie.exceptions import MissingAliasError

        calls = {}

        def fake_pull(*args, **kwargs):
            calls["args"] = args
            calls["kwargs"] = kwargs
            return None

        monkeypatch.setattr(refgenie_minimal, "pull", fake_pull)

        with pytest.raises(MissingAliasError):
            refgenie_minimal.build_asset(
                recipe_name="fasta",
                genome_name="nonexistent_genome",
                asset_group_name="fasta",
                pull_parents=True,
            )

        assert calls["args"] == ()
        assert calls["kwargs"] == {
            "asset_group_name": "fasta",
            "alias_name": "nonexistent_genome",
        }


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
    def test_estimate_pull_size(self, refgenie_minimal, asset_list, expected):
        assert refgenie_minimal.asset.estimate_pull_size(asset_list) == expected


class TestPullMultiple:
    """pull_multiple continues past individual failures and skips invalid entries."""

    def test_continues_on_individual_failure(self, refgenie_minimal):
        asset_list = [
            {
                "genome_digest": "abc123",
                "asset_group_name": "fasta",
                "asset_name": "default",
                "server_url": "http://example.com",
            },
            {
                "genome_digest": "def456",
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

        with patch.object(refgenie_minimal.asset._asset_puller, "pull", side_effect=mock_pull):
            results = refgenie_minimal.asset.pull_multiple(asset_list)
        assert len(results) == 1

    def test_skips_entries_missing_required_fields(self, refgenie_minimal):
        asset_list = [
            {"genome_digest": None, "asset_group_name": "fasta"},
            {"genome_digest": "abc123"},  # missing asset_group_name
        ]
        assert refgenie_minimal.asset.pull_multiple(asset_list) == []


# ===========================================================================
# list-remote: querying subscribed servers for available assets
# ===========================================================================


class TestListRemote:
    """Test list-remote (querying remote servers for available assets).

    Component tier: these build real genome folders, asset files and archives
    on disk alongside SQLite. Deselected from the bare `pytest` inner loop;
    run with `pytest -m component`.
    """

    pytestmark = pytest.mark.component

    @pytest.fixture
    def list_remote_client(self, engine, tmp_path, fixtures_path):
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
        client_rg.configuration.subscribe(mock_server_url)

        with serve_refgenie(client_rg, server_rg, mock_server_url):
            yield client_rg, mock_server_url

    def test_list_remote_returns_assets_and_aliases(self, list_remote_client):
        """list-remote returns both asset and alias data from a subscribed server."""
        client_rg, mock_server_url = list_remote_client
        asset_data, aliases_data = client_rg.asset.list_remote()

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

    def test_assets_table_remote(self, list_remote_client):
        """asset.remote_table() returns Rich tables."""
        from rich.table import Table

        client_rg, _ = list_remote_client
        tables = client_rg.asset.remote_table()

        assert isinstance(tables, list)
        assert len(tables) > 0
        for table in tables:
            assert isinstance(table, Table)


def test_list_remote_no_subscriptions_returns_empty(refgenie_minimal):
    """list-remote with no subscriptions returns empty dicts."""
    asset_data, aliases_data = refgenie_minimal.asset.list_remote()
    assert asset_data == {}
    assert aliases_data == {}


class TestOverwriteConfirmation:
    """The "replace the existing asset directory?" prompt goes through the
    caller's confirmer, not a bare `rich.prompt.Confirm.ask`.

    That branch used to call rich directly, which meant `confirm=deny` did not
    guard it -- only `force=False` happened to avoid reaching it. A pull driven
    from a server request or a worker thread would block there forever, on a
    terminal nobody is watching.
    """

    @pytest.fixture
    def puller_at_the_prompt(self, mock_puller):
        """A puller positioned exactly on the "directory already exists" branch."""
        from refgenie.managers.asset.puller import PulledAssetMetadata

        asset_digest = "existingassetdigest"
        mock_puller._fetch_asset_metadata = MagicMock(
            return_value=PulledAssetMetadata(
                asset_metadata={},
                asset_group_metadata={},
                asset_digest=asset_digest,
                declared_modes=["archive"],
                asset_parents=[],
                asset_name="default",
                asset_class_name="fasta",
                asset_names=[],
            )
        )
        mock_puller.get_client = MagicMock(return_value=MagicMock())
        asset_dir = mock_puller.data_folder / "abc123digest" / "fasta" / asset_digest
        asset_dir.mkdir(parents=True)
        existing_row = MagicMock()
        existing_row.path = str(asset_dir)
        mock_puller._asset_manager.get_by_digest.return_value = existing_row
        return mock_puller

    def test_the_callers_confirmer_decides(self, puller_at_the_prompt):
        asked = []

        def refuse(message):
            asked.append(message)
            return False

        with patch("refgenie.utils.prompt.ask", side_effect=AssertionError("read stdin")):
            result = puller_at_the_prompt.pull(
                asset_group_name="fasta", alias_name="rCRSd", force=None, confirm=refuse
            )

        assert result is None, "a declined overwrite is a skip, not a pull"
        assert asked and "Replace existing" in asked[0]
