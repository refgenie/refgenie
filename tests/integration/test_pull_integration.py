"""
Integration tests for pulling across multiple remotes and serving modes.

Every test here talks to a real `refgenie serve` subprocess over a real socket:

- Multi-remote pull (rCRSd from server A, demo from server B)
- File-mode pull (asset-level serving_modes_override forces per-file downloads)
- Cross-mode content fidelity (file-served bytes == archive-served bytes)

Route-level coverage of the server API belongs in the component tier
(`tests/test_server.py`); this module exists only for what needs a
real process. There are no in-process apps and no hand-rolled route fakes here.

All servers use SQLite databases (no PostgreSQL dependency).

Run via: ./tests/scripts/test-integration.sh
"""

# Pull tests need separate subprocess servers (they require a real `refgenie
# serve` process so clients can fetch its OpenAPI spec). These are inherently
# expensive; minimize the number of them.
#
# Use `build_server_env()` (tests/integration/conftest.py) to add multiple
# genomes with different serving_modes_override values to a single DB, rather
# than creating separate DBs per mode.
#
# See tests/README.md for general integration-test fixture conventions.

import os
import tarfile
from pathlib import Path

import httpx
import pytest

from tests.integration.conftest import (
    _create_client_env,
    normalize_fasta,
    run_refgenie,
    RCRSD_FASTA,
)

pytestmark = pytest.mark.skipif(
    os.getenv("RUN_INTEGRATION_TESTS") != "true",
    reason="Integration tests disabled. Run ./tests/scripts/test-integration.sh",
)


@pytest.mark.shared_state
class TestMultiRemotePull:
    """Pull from two servers: rCRSd from server A, demo from server B."""

    def test_pull_from_different_servers(self, multi_remote_servers, tmp_path):
        """Client subscribed to two servers can pull different genomes from each."""
        server_a = multi_remote_servers["server_a"]
        server_b = multi_remote_servers["server_b"]

        client_env = _create_client_env(tmp_path, server_a["url"], server_b["url"])

        # Pull rCRSd (should come from server A)
        result_a = run_refgenie(
            "pull", "-g", "rCRSd", "--all", "--force", env=client_env, check=False
        )
        assert result_a.returncode == 0, (
            f"Pull rCRSd failed:\nSTDOUT: {result_a.stdout[-500:]}\n"
            f"STDERR: {result_a.stderr[-500:]}"
        )

        # Pull demo (should come from server B)
        result_b = run_refgenie(
            "pull", "-g", "demo", "--all", "--force", env=client_env, check=False
        )
        assert result_b.returncode == 0, (
            f"Pull demo failed:\nSTDOUT: {result_b.stdout[-500:]}\nSTDERR: {result_b.stderr[-500:]}"
        )

        # Verify both assets exist via seek
        for genome in ["rCRSd", "demo"]:
            seek_result = run_refgenie(
                "seek", f"{genome}/fasta:default", env=client_env, check=False
            )
            assert seek_result.returncode == 0, (
                f"Seek {genome} failed:\nSTDOUT: {seek_result.stdout[-500:]}\n"
                f"STDERR: {seek_result.stderr[-500:]}"
            )
            assert os.path.isfile(seek_result.stdout.strip()), f"Asset not found for {genome}"


@pytest.mark.shared_state
class TestFileModePull:
    """Asset-level serving_modes_override reaches staging, and the resulting
    file-only server is pulled from by downloading individual files."""

    def test_pull_file_mode_asset(self, file_mode_pull_server, tmp_path):
        """File-mode pull downloads individual files and produces a valid asset."""
        asset_digest = file_mode_pull_server["genomes"]["rCRSd"]["asset_digest"]
        url = file_mode_pull_server["url"]

        resp = httpx.get(f"{url}/v4/assets/{asset_digest}")
        assert resp.status_code == 200, resp.text
        # Asset-level serving_modes_override reached staging and is surfaced to the
        # puller; without this the client would try the archive path.
        assert resp.json()["serving_modes"] == ["file"]
        assert httpx.get(f"{url}/v4/archives/{asset_digest}/download").status_code == 404

        client_env = _create_client_env(tmp_path, url)

        result = run_refgenie(
            "pull", "-g", "rCRSd", "--all", "--force", env=client_env, check=False
        )
        assert result.returncode == 0, (
            f"File-mode pull failed:\nSTDOUT: {result.stdout[-1000:]}\n"
            f"STDERR: {result.stderr[-1000:]}"
        )

        # Verify asset exists
        seek_result = run_refgenie("seek", "rCRSd/fasta:default", env=client_env, check=False)
        assert seek_result.returncode == 0, (
            f"Seek failed:\nPull STDOUT: {result.stdout[-500:]}\n"
            f"Seek STDERR: {seek_result.stderr[-500:]}"
        )
        seek_path = seek_result.stdout.strip()
        assert os.path.isfile(seek_path), f"Pulled asset not found at {seek_path}"

        # Verify content matches original
        pulled = normalize_fasta(Path(seek_path).read_text())
        original = normalize_fasta(RCRSD_FASTA.read_text())
        assert pulled == original


@pytest.mark.shared_state
class TestCrossModeContentFidelity:
    """The one assertion nothing else in the suite makes: for a both-modes asset,
    the bytes served by the file endpoint are identical to the bytes inside the
    tarball served by the archive endpoint."""

    def test_content_consistent_between_modes(self, multi_remote_servers, tmp_path):
        """FASTA from the file endpoint matches FASTA extracted from the archive."""
        server_b = multi_remote_servers["server_b"]
        asset_digest = server_b["genomes"]["demo"]["asset_digest"]
        url = server_b["url"]

        # File endpoint: {"asset_digest": ..., "files": [...]}, then raw bytes.
        listing = httpx.get(f"{url}/v4/assets/{asset_digest}/files")
        assert listing.status_code == 200, listing.text
        files = listing.json()["files"]
        fa_file = next(f for f in files if f.endswith(".fa"))
        file_resp = httpx.get(f"{url}/v4/assets/{asset_digest}/files/{fa_file}")
        assert file_resp.status_code == 200, file_resp.text
        file_content = file_resp.content

        # Archive endpoint: no http/https remote is configured, so the real route
        # serves the tarball locally rather than redirecting.
        archive_resp = httpx.get(f"{url}/v4/archives/{asset_digest}/download")
        assert archive_resp.status_code == 200, archive_resp.text
        assert archive_resp.content[:2] == b"\x1f\x8b", "Expected gzip magic bytes"

        archive_file = tmp_path / "asset.tgz"
        archive_file.write_bytes(archive_resp.content)
        with tarfile.open(archive_file, "r:gz") as tar:
            tar.extractall(path=tmp_path / "extracted")

        fasta_files = list((tmp_path / "extracted").rglob("*.fa"))
        assert len(fasta_files) > 0, "No FASTA file in the extracted archive"
        archive_content = fasta_files[0].read_bytes()

        assert file_content == archive_content
