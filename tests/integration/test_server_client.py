"""
Cross-package integration tests for refgenie + refgenconf + refgenieserver.

Tests the full client->server->config chain using the REAL refgenieserver:
- Server endpoints respond correctly (exercising actual refgenieserver code)
- RefGenConf.listr() discovers remote assets
- RefGenConf.pull() downloads and unpacks assets
- Alias operations work end-to-end
- Subscribe/unsubscribe operations work
- Config survives write/reload cycles
- Known bugs are caught by regression tests
"""

from __future__ import annotations

import hashlib
import io
import os
import subprocess
import sys
import tarfile

import pytest
import requests

from .conftest import GENOME_ALIAS, GENOME_DESC, GENOME_DIGEST, FASTA_CONTENT


# =============================================================================
# Test group 1: Server endpoint tests (via direct HTTP)
# =============================================================================


class TestServerEndpoints:
    """Verify the real refgenieserver serves assets correctly via HTTP."""

    def test_genomes_list(self, test_server_url):
        """Server lists available genomes."""
        resp = requests.get(f"{test_server_url}/v3/genomes/list")
        assert resp.status_code == 200
        data = resp.json()
        assert isinstance(data, list)
        assert GENOME_DIGEST in data

    def test_alias_dict(self, test_server_url):
        """Server returns alias-to-digest mapping."""
        resp = requests.get(f"{test_server_url}/v3/genomes/alias_dict")
        assert resp.status_code == 200
        data = resp.json()
        assert GENOME_DIGEST in data
        assert GENOME_ALIAS in data[GENOME_DIGEST]

    def test_assets_list(self, test_server_url):
        """Server lists assets for genomes."""
        resp = requests.get(f"{test_server_url}/v3/assets/list")
        assert resp.status_code == 200
        data = resp.json()
        assert GENOME_DIGEST in data
        assert "fasta" in data[GENOME_DIGEST]

    def test_assets_list_with_seek_keys(self, test_server_url):
        """Server lists assets with seek keys when requested."""
        resp = requests.get(
            f"{test_server_url}/v3/assets/list", params={"includeSeekKeys": True}
        )
        assert resp.status_code == 200
        data = resp.json()
        assert GENOME_DIGEST in data
        # With includeSeekKeys=True, the real server returns tags in the list
        genome_data = data[GENOME_DIGEST]
        # Should contain fasta-related entries
        if isinstance(genome_data, list):
            assert any("fasta" in item for item in genome_data)
        else:
            assert "fasta" in genome_data

    def test_asset_attributes(self, test_server_url):
        """Server returns asset attributes."""
        resp = requests.get(f"{test_server_url}/v3/assets/attrs/{GENOME_DIGEST}/fasta")
        assert resp.status_code == 200
        data = resp.json()
        assert "asset_digest" in data
        assert "archive_digest" in data
        assert "archive_size" in data
        assert "seek_keys" in data

    def test_genome_attributes(self, test_server_url):
        """Server returns genome attributes."""
        resp = requests.get(f"{test_server_url}/v3/genomes/attrs/{GENOME_DIGEST}")
        assert resp.status_code == 200
        data = resp.json()
        assert "genome_description" in data
        assert data["genome_description"] == GENOME_DESC

    def test_default_tag(self, test_server_url):
        """Server returns default tag for an asset."""
        resp = requests.get(
            f"{test_server_url}/v3/assets/default_tag/{GENOME_DIGEST}/fasta"
        )
        assert resp.status_code == 200
        assert resp.text == "default"

    def test_asset_digest(self, test_server_url):
        """Server returns asset digest."""
        resp = requests.get(
            f"{test_server_url}/v3/assets/asset_digest/{GENOME_DIGEST}/fasta"
        )
        assert resp.status_code == 200
        assert len(resp.text) == 32  # md5 hex digest

    def test_archive_download(self, test_server_url):
        """Server serves asset archive for download."""
        resp = requests.get(
            f"{test_server_url}/v3/assets/archive/{GENOME_DIGEST}/fasta",
            params={"tag": "default"},
        )
        assert resp.status_code == 200
        # Verify it's a valid gzip file (starts with magic bytes)
        assert resp.content[:2] == b"\x1f\x8b"

    def test_genome_alias_resolution(self, test_server_url):
        """Server resolves genome alias to digest."""
        resp = requests.get(
            f"{test_server_url}/v3/genomes/genome_digest/{GENOME_ALIAS}"
        )
        assert resp.status_code == 200
        assert resp.text == GENOME_DIGEST

    def test_openapi_spec(self, test_server_url):
        """Server provides valid OpenAPI spec with correct operation IDs."""
        resp = requests.get(f"{test_server_url}/openapi.json")
        assert resp.status_code == 200
        spec = resp.json()
        assert "openapi" in spec
        assert "paths" in spec

        # Verify key operation IDs are present (from real refgenieserver)
        operation_ids = set()
        for path_data in spec["paths"].values():
            for method_data in path_data.values():
                if "operationId" in method_data:
                    operation_ids.add(method_data["operationId"])

        assert "v3custom_Id_assets" in operation_ids
        assert "v3custom_Id_archive" in operation_ids

    def test_safely_get_example_returns_real_values(self, server_rgc):
        """safely_get_example returns actual values from config, not fallback defaults."""
        from refgenieserver.helpers import safely_get_example

        # This should return a real genome from the config, not the fallback
        result = safely_get_example(
            server_rgc,
            "genome",
            "genomes_list",
            "fallback_default_value",
        )
        # If __getattr__ bug exists, this returns the fallback
        assert result != "fallback_default_value", (
            "safely_get_example fell back to default - likely __getattr__ bug"
        )
        # Should return the actual genome alias (genomes_list returns aliases)
        assert result == GENOME_ALIAS


# =============================================================================
# Test group 2: RefGenConf client operations against test server
# =============================================================================


class TestRefGenConfClientOperations:
    """Test RefGenConf client methods against the real refgenieserver."""

    def test_list_remote(self, client_rgc):
        """RefGenConf.listr() returns assets from test server."""
        remote_data = client_rgc.listr()
        assert len(remote_data) > 0
        # listr returns {server_url: {genome_alias: assets}}
        for server_url, genome_data in remote_data.items():
            assert len(genome_data) > 0
            # Should contain the genome alias
            assert GENOME_ALIAS in genome_data

    def test_list_remote_as_digests(self, client_rgc):
        """RefGenConf.listr(as_digests=True) returns digests instead of aliases."""
        remote_data = client_rgc.listr(as_digests=True)
        assert len(remote_data) > 0
        for server_url, genome_data in remote_data.items():
            assert GENOME_DIGEST in genome_data

    def test_list_remote_filtered_by_genome(self, client_rgc):
        """RefGenConf.listr(genome=...) filters by genome."""
        remote_data = client_rgc.listr(genome=GENOME_ALIAS)
        assert len(remote_data) > 0
        for server_url, genome_data in remote_data.items():
            assert GENOME_ALIAS in genome_data

    def test_pull_asset(self, client_rgc):
        """RefGenConf.pull() downloads and unpacks an asset from the test server."""
        result = client_rgc.pull(GENOME_ALIAS, "fasta", "default", force=True)
        assert result is not None
        gat, archive_data, server_url = result
        assert gat[1] == "fasta"
        assert gat[2] == "default"
        assert archive_data is not None
        assert server_url is not None

    def test_pull_creates_files(self, client_rgc):
        """After pull, asset files exist on disk."""
        client_rgc.pull(GENOME_ALIAS, "fasta", "default", force=True)
        seek_path = client_rgc.seek(GENOME_ALIAS, "fasta")
        assert os.path.exists(seek_path)

    def test_pull_updates_config(self, client_rgc):
        """After pull, config reflects the new asset."""
        client_rgc.pull(GENOME_ALIAS, "fasta", "default", force=True)
        local_assets = client_rgc.list()
        # The genome (by digest or alias) should now appear in local assets
        assert len(local_assets) > 0

    def test_seek_after_pull(self, client_rgc):
        """seek() returns valid path after pulling an asset."""
        client_rgc.pull(GENOME_ALIAS, "fasta", "default", force=True)
        path = client_rgc.seek(GENOME_ALIAS, "fasta")
        assert os.path.isfile(path)
        # Verify the file has content
        assert os.path.getsize(path) > 0

    def test_seek_fasta_seek_key_after_pull(self, client_rgc):
        """seek() with seek_key='fasta' returns the FASTA file path."""
        client_rgc.pull(GENOME_ALIAS, "fasta", "default", force=True)
        path = client_rgc.seek(GENOME_ALIAS, "fasta", seek_key="fasta")
        assert os.path.isfile(path)
        assert path.endswith(".fa")

    def test_seek_dir_after_pull(self, client_rgc):
        """seek() with seek_key='dir' returns the asset directory."""
        client_rgc.pull(GENOME_ALIAS, "fasta", "default", force=True)
        path = client_rgc.seek(GENOME_ALIAS, "fasta", seek_key="dir")
        assert os.path.isdir(path)


# =============================================================================
# Test group 3: Alias operations
# =============================================================================


class TestAliasOperations:
    """Test genome alias operations via the client."""

    def test_genome_alias_digest(self, client_rgc):
        """After pull, get_genome_alias_digest resolves alias to digest."""
        client_rgc.pull(GENOME_ALIAS, "fasta", "default", force=True)
        digest = client_rgc.get_genome_alias_digest(alias=GENOME_ALIAS)
        assert digest == GENOME_DIGEST

    def test_get_genome_alias(self, client_rgc):
        """After pull, get_genome_alias returns alias for digest."""
        client_rgc.pull(GENOME_ALIAS, "fasta", "default", force=True)
        alias = client_rgc.get_genome_alias(digest=GENOME_DIGEST)
        # Should be the genome alias or a list containing it
        if isinstance(alias, list):
            assert GENOME_ALIAS in alias
        else:
            assert alias == GENOME_ALIAS

    def test_set_genome_alias(self, client_rgc):
        """Set a custom genome alias and verify it resolves."""
        client_rgc.pull(GENOME_ALIAS, "fasta", "default", force=True)
        client_rgc.set_genome_alias(
            genome="my_mito", digest=GENOME_DIGEST, no_write=True
        )
        # Verify alias resolves to the correct digest
        digest = client_rgc.get_genome_alias_digest("my_mito")
        assert digest == GENOME_DIGEST

    def test_genome_aliases_property(self, client_rgc):
        """genome_aliases property returns alias mapping."""
        client_rgc.pull(GENOME_ALIAS, "fasta", "default", force=True)
        aliases = client_rgc.genome_aliases
        assert isinstance(aliases, dict)
        assert len(aliases) > 0


# =============================================================================
# Test group 4: Subscribe/unsubscribe
# =============================================================================


class TestSubscription:
    """Test server subscription management."""

    def test_subscribe(self, client_rgc):
        """Subscribe adds server to config."""
        from refgenconf.const import CFG_SERVERS_KEY

        client_rgc.subscribe(["http://test-server.example.com"], no_write=True)
        assert "http://test-server.example.com" in client_rgc[CFG_SERVERS_KEY]

    def test_unsubscribe(self, client_rgc):
        """Unsubscribe removes server from config."""
        from refgenconf.const import CFG_SERVERS_KEY

        client_rgc.subscribe(["http://test-server.example.com"], no_write=True)
        client_rgc.unsubscribe(["http://test-server.example.com"], no_write=True)
        assert "http://test-server.example.com" not in client_rgc[CFG_SERVERS_KEY]

    def test_subscribe_deduplicates(self, client_rgc):
        """Subscribing to the same server twice doesn't duplicate."""
        from refgenconf.const import CFG_SERVERS_KEY

        original_count = len(client_rgc[CFG_SERVERS_KEY])
        server = client_rgc[CFG_SERVERS_KEY][0]
        client_rgc.subscribe([server], no_write=True)
        assert len(client_rgc[CFG_SERVERS_KEY]) == original_count


# =============================================================================
# Test group 5: Config persistence
# =============================================================================


class TestConfigPersistence:
    """Test that config changes survive write/reload cycles."""

    def test_write_and_reload(self, client_rgc, tmp_path):
        """Config survives write/reload cycle after pull."""
        import yaml as _yaml

        from refgenconf import RefGenConf

        client_rgc.pull(GENOME_ALIAS, "fasta", "default", force=True)

        # Export config to YAML and write to a new path (avoids lock issues)
        new_cfg_path = str(tmp_path / "saved_config.yaml")
        with open(new_cfg_path, "w") as f:
            f.write(client_rgc.to_yaml())

        # Reload from disk
        reloaded = RefGenConf.from_yaml_file(new_cfg_path)
        original_path = client_rgc.seek(GENOME_ALIAS, "fasta")
        reloaded_path = reloaded.seek(GENOME_ALIAS, "fasta")
        assert reloaded_path == original_path

    def test_pull_persists_to_disk(self, client_rgc):
        """Pull writes config changes to disk automatically."""
        from refgenconf import RefGenConf

        client_rgc.pull(GENOME_ALIAS, "fasta", "default", force=True)

        # Reload from the same config file
        if client_rgc.file_path:
            reloaded = RefGenConf.from_yaml_file(client_rgc.file_path)
            assets = reloaded.list()
            assert len(assets) > 0


# =============================================================================
# Test group 6: Refgenie CLI integration
# =============================================================================


class TestRefgenieCLI:
    """Test refgenie CLI commands against the integration setup."""

    def test_list_command(self, client_config_path):
        """refgenie list shows config info."""
        result = subprocess.run(
            [sys.executable, "-m", "refgenie", "list", "-c", client_config_path],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0

    def test_listr_command(self, client_config_path):
        """refgenie listr lists remote assets."""
        result = subprocess.run(
            [sys.executable, "-m", "refgenie", "listr", "-c", client_config_path],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        # Output should mention genome alias or the genome digest
        assert GENOME_ALIAS in result.stdout or GENOME_DIGEST in result.stdout

    def test_pull_command(self, client_config_path):
        """refgenie pull downloads an asset via CLI."""
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "refgenie",
                "pull",
                "-c",
                client_config_path,
                f"{GENOME_ALIAS}/fasta:default",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0

    def test_seek_after_pull_command(self, client_config_path):
        """refgenie seek returns path after pull."""
        # First pull
        subprocess.run(
            [
                sys.executable,
                "-m",
                "refgenie",
                "pull",
                "-c",
                client_config_path,
                f"{GENOME_ALIAS}/fasta:default",
            ],
            capture_output=True,
            text=True,
        )
        # Then seek
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "refgenie",
                "seek",
                "-c",
                client_config_path,
                f"{GENOME_ALIAS}/fasta",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        # Output should be a file path
        path = result.stdout.strip()
        assert len(path) > 0


# =============================================================================
# Test group 7: Bug regression tests
# =============================================================================


class TestKnownBugRegressions:
    """Tests that catch specific bugs we know existed."""

    def test_server_startup_no_context_manager_error(self, test_server_url):
        """The server starts without NotImplementedError from 'with rgc as r:'.

        The old bug was that main.py:61 used 'with rgc as r:' which
        raises NotImplementedError in yacman v1's YAMLConfigManager.
        If we got here, the server started successfully.
        """
        resp = requests.get(f"{test_server_url}/v3/genomes/list")
        assert resp.status_code == 200

    def test_refgenie_init_no_skip_read_lock_error(self, tmp_path):
        """refgenie init does not crash with unexpected skip_read_lock param.

        The old bug was that cli.py:135 passed skip_read_lock=... to
        RefGenConf() which no longer accepts that parameter.
        """
        config_path = str(tmp_path / "test_init.yaml")
        result = subprocess.run(
            [sys.executable, "-m", "refgenie", "init", "-c", config_path],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, f"stderr: {result.stderr}"


# =============================================================================
# Test group 8: Smoke test - verify real refgenieserver code
# =============================================================================


class TestRealServerVerification:
    """Verify the test server is the real refgenieserver, not a mock."""

    def test_server_is_real_refgenieserver(self, test_server_url):
        """Verify the test server is the real refgenieserver, not a mock."""
        resp = requests.get(f"{test_server_url}/openapi.json")
        spec = resp.json()
        # Real refgenieserver uses PKG_NAME = "refgenieserver" as the title
        assert spec["info"]["title"] == "refgenieserver"
        # Real refgenieserver has specific API version tags
        tag_names = [t["name"] for t in spec.get("tags", [])]
        assert "APIv3" in tag_names

    def test_private_api_genomes_dict(self, test_server_url):
        """The private genomes/dict endpoint returns full genome config."""
        resp = requests.get(f"{test_server_url}/_private_api/genomes/dict")
        assert resp.status_code == 200
        data = resp.json()
        assert GENOME_DIGEST in data
        genome_data = data[GENOME_DIGEST]
        assert "assets" in genome_data
        assert "fasta" in genome_data["assets"]

    def test_private_api_tag_fields_complete(self, test_server_url):
        """The private API returns all fields required by the Tag data model."""
        resp = requests.get(f"{test_server_url}/_private_api/genomes/dict")
        data = resp.json()
        tag = data[GENOME_DIGEST]["assets"]["fasta"]["tags"]["default"]
        required_fields = [
            "asset_path",
            "asset_digest",
            "archive_digest",
            "asset_size",
            "archive_size",
            "seek_keys",
            "asset_parents",
            "asset_children",
        ]
        for field in required_fields:
            assert field in tag, f"Missing required Tag field: {field}"


# =============================================================================
# Test group 9: End-to-end pull verification
# =============================================================================


class TestEndToEndPullVerification:
    """Verify that pulled assets contain the correct data."""

    def test_pull_fasta_content_matches(self, client_rgc):
        """After pull, the FASTA file content matches the original test data."""
        client_rgc.pull(GENOME_ALIAS, "fasta", "default", force=True)
        fasta_path = client_rgc.seek(GENOME_ALIAS, "fasta")
        with open(fasta_path) as f:
            content = f.read()
        # The test data is t7.fa.gz (bacteriophage T7, accession V01146.1)
        assert content == FASTA_CONTENT

    def test_pull_fasta_file_not_empty(self, client_rgc):
        """After pull, the FASTA file is non-empty and starts with a header."""
        client_rgc.pull(GENOME_ALIAS, "fasta", "default", force=True)
        fasta_path = client_rgc.seek(GENOME_ALIAS, "fasta")
        with open(fasta_path) as f:
            content = f.read()
        assert len(content) > 0
        assert content.startswith(">")

    def test_pull_fasta_seek_key_matches_default(self, client_rgc):
        """seek() and seek(seek_key='fasta') point to the same file."""
        client_rgc.pull(GENOME_ALIAS, "fasta", "default", force=True)
        default_path = client_rgc.seek(GENOME_ALIAS, "fasta")
        fasta_path = client_rgc.seek(GENOME_ALIAS, "fasta", seek_key="fasta")
        # Both should resolve to the same actual file
        assert os.path.samefile(default_path, fasta_path)

    def test_cli_pull_fasta_content_matches(self, client_config_path):
        """After CLI pull, the FASTA content matches the original test data."""
        subprocess.run(
            [
                sys.executable,
                "-m",
                "refgenie",
                "pull",
                "-c",
                client_config_path,
                f"{GENOME_ALIAS}/fasta:default",
            ],
            capture_output=True,
            text=True,
            check=True,
        )
        result = subprocess.run(
            [
                sys.executable,
                "-m",
                "refgenie",
                "seek",
                "-c",
                client_config_path,
                f"{GENOME_ALIAS}/fasta",
            ],
            capture_output=True,
            text=True,
            check=True,
        )
        fasta_path = result.stdout.strip()
        with open(fasta_path) as f:
            content = f.read()
        # The test data is t7.fa.gz (bacteriophage T7, accession V01146.1)
        assert content == FASTA_CONTENT


# =============================================================================
# Test group 10: Archive integrity
# =============================================================================


class TestArchiveIntegrity:
    """Verify the archive format is correct."""

    def test_archive_is_valid_tgz(self, test_server_url):
        """The archive endpoint returns valid gzip data."""
        resp = requests.get(
            f"{test_server_url}/v3/assets/archive/{GENOME_DIGEST}/fasta",
            params={"tag": "default"},
        )
        assert resp.status_code == 200
        assert resp.content[:2] == b"\x1f\x8b"  # gzip magic bytes

    def test_archive_extracts_correctly(self, test_server_url, tmp_path):
        """The archive can be extracted and contains the expected files."""
        resp = requests.get(
            f"{test_server_url}/v3/assets/archive/{GENOME_DIGEST}/fasta",
            params={"tag": "default"},
        )
        with tarfile.open(fileobj=io.BytesIO(resp.content), mode="r:gz") as tar:
            names = tar.getnames()
        # Should contain the asset files under "default/" directory
        assert any("default" in n for n in names)

    def test_archive_checksum_matches(self, test_server_url):
        """The archive's computed checksum matches what the server reports."""
        # Get reported digest
        resp_digest = requests.get(
            f"{test_server_url}/v3/assets/archive_digest/{GENOME_DIGEST}/fasta"
        )
        reported_digest = resp_digest.text
        # Download archive and compute digest
        resp_archive = requests.get(
            f"{test_server_url}/v3/assets/archive/{GENOME_DIGEST}/fasta",
            params={"tag": "default"},
        )
        computed_digest = hashlib.md5(resp_archive.content).hexdigest()
        assert computed_digest == reported_digest
