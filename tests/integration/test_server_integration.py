"""
Server Integration Tests -- HTTP/server tests against a real PostgreSQL.

This file covers server endpoint and data channel scenarios.

Run via: ./tests/scripts/test-integration.sh
"""

import os
import urllib.request
from pathlib import Path

import httpx
import pytest
import yaml

from tests.helpers import TESTS_DATA_DIR, register_fasta


# Skip all tests in this module unless RUN_INTEGRATION_TESTS=true
pytestmark = pytest.mark.skipif(
    os.getenv("RUN_INTEGRATION_TESTS") != "true",
    reason="Integration tests disabled. Run ./tests/scripts/test-integration.sh",
)


# =============================================================================
# Scenario 11: Server Endpoints
# =============================================================================


class TestServerEndpoints:
    """The real app boots and answers against PostgreSQL.

    Every endpoint's contract is asserted at the unit/component tier
    (tests/test_server.py). This tier only proves the app comes
    up against a real PostgreSQL engine rather than in-memory SQLite.
    """

    def test_app_boots_against_postgres(self, client):
        assert client.get("/v4/healthcheck").json()["status"] == "ok"
        assert "genomes" in client.get("/v4/summary").json()


# =============================================================================
# Scenario 12: Data Channel Integration
# =============================================================================


class TestDataChannelIntegration:
    """Scenario 12: HTTP data channel server integration."""

    def test_data_channel_integration(self, refgenie_instance, data_channel_url, fixtures_path):
        """HTTP data channel: server running, serves files, add/remove channel."""
        # -- Server accessible --
        try:
            with urllib.request.urlopen(f"{data_channel_url}/", timeout=5) as response:
                assert response.status == 200
        except Exception as e:
            pytest.fail(f"Data channel server not accessible at {data_channel_url}: {e}")

        # -- Serves asset class YAML --
        if (fixtures_path / "fasta_asset_class.yaml").exists():
            url = f"{data_channel_url}/fasta_asset_class.yaml"
            with urllib.request.urlopen(url, timeout=5) as response:
                content = response.read().decode("utf-8")
                assert "fasta" in content.lower()
                assert response.status == 200

        # -- Serves recipe YAML --
        if (fixtures_path / "fasta_asset_recipe.yaml").exists():
            url = f"{data_channel_url}/fasta_asset_recipe.yaml"
            with urllib.request.urlopen(url, timeout=5) as response:
                content = response.read().decode("utf-8")
                assert "fasta" in content.lower()
                assert response.status == 200

        # -- Add channel from HTTP --
        from refgenie.db.tables import DataChannelType

        channel_name = "test_integration_channel"

        # Remove if exists from previous run
        refgenie_instance.sources.remove_channel(channel_name)

        # Add channel
        refgenie_instance.sources.add_channel(
            name=channel_name,
            type=DataChannelType.http,
            index_address=f"{data_channel_url}/",
            description="Integration test HTTP channel",
        )

        # Verify added
        channels = list(refgenie_instance.sources.list_channels())
        channel_names = [c.name for c in channels]
        assert channel_name in channel_names

        # Cleanup
        refgenie_instance.sources.remove_channel(channel_name)


# TestDownloadWithProgress (download progress bar handling) lives in
# tests/test_server_client.py -- it tests client behavior with an in-process
# TestClient and needs neither Docker nor postgres.


# --- merged from test_store_backed_serving.py --------------------------------
# Integration tests for store-backed server mode: refgenie's production serving
# path, a remote RefgetStore served over HTTP with the refgenie server in
# store-backed mode (REFGENIE_REFGET_STORE_URL set). The `store_backed_server`
# fixture (in conftest.py) builds a small on-disk RefgetStore from refget test
# FASTAs, serves it over HTTP via a `python -m http.server` subprocess, and
# starts `refgenie serve` as a subprocess pointed at the HTTP store.


@pytest.mark.shared_state
class TestStoreBackedServing:
    """The refgenie server in store-backed mode serves store-derived seqcol data."""

    def test_service_info_reports_store_enabled(self, store_backed_server):
        """GET /seqcol/service-info reflects the store-backed backend and counts."""
        resp = httpx.get(f"{store_backed_server['url']}/seqcol/service-info", timeout=30)
        assert resp.status_code == 200
        data = resp.json()
        assert data["id"] == "org.refgenie.seqcol"
        store_info = data["seqcol"]["refget_store"]
        assert store_info["enabled"] is True
        assert store_info["url"] == store_backed_server["store_url"]
        # 3 collections were loaded (base.fa, subset.fa, different_names.fa).
        assert store_info["n_collections"] == 3
        # All three collections reuse a shared pool of 3 distinct sequences.
        assert store_info["n_sequences"] == 3
        assert store_info["backend_type"] == "refget_store"

    def test_list_collections(self, store_backed_server):
        """GET /seqcol/list/collection lists the store's collection digests."""
        resp = httpx.get(f"{store_backed_server['url']}/seqcol/list/collection", timeout=30)
        assert resp.status_code == 200
        data = resp.json()
        assert data["pagination"]["total"] == 3
        results = set(data["results"])
        assert results == set(store_backed_server["digests"])

    def test_get_collection(self, store_backed_server):
        """GET /seqcol/collection/{digest} returns a valid seqcol document."""
        digest = store_backed_server["digests"][0]
        resp = httpx.get(f"{store_backed_server['url']}/seqcol/collection/{digest}", timeout=30)
        assert resp.status_code == 200
        data = resp.json()
        # Core seqcol attributes derived from the store.
        for key in ("names", "lengths", "sequences"):
            assert key in data, f"missing {key} in collection document"
        assert len(data["names"]) == len(data["lengths"]) == len(data["sequences"])
        assert len(data["names"]) > 0
        assert all(s.startswith("SQ.") for s in data["sequences"])

    def test_get_collection_level2(self, store_backed_server):
        """GET /seqcol/collection/{digest}?level=2 matches the default level."""
        digest = store_backed_server["digests"][0]
        url = store_backed_server["url"]
        default = httpx.get(f"{url}/seqcol/collection/{digest}", timeout=30).json()
        level2 = httpx.get(f"{url}/seqcol/collection/{digest}?level=2", timeout=30)
        assert level2.status_code == 200
        assert level2.json() == default

    def test_get_collection_unknown_digest_404(self, store_backed_server):
        """An unknown collection digest returns 404."""
        url = store_backed_server["url"]
        resp = httpx.get(f"{url}/seqcol/collection/NONEXISTENT_DIGEST_12345678901234", timeout=30)
        assert resp.status_code == 404

    def test_comparison(self, store_backed_server):
        """GET /seqcol/comparison/{a}/{b} compares two store collections."""
        digests = store_backed_server["digests"]
        a, b = digests[0], digests[1]
        resp = httpx.get(f"{store_backed_server['url']}/seqcol/comparison/{a}/{b}", timeout=30)
        assert resp.status_code == 200
        data = resp.json()
        assert data["digests"] == {"a": a, "b": b}
        assert "attributes" in data
        assert "array_elements" in data

    def test_drs_service_info(self, store_backed_server):
        """GET /ga4gh/drs/service-info responds in store-backed mode."""
        resp = httpx.get(f"{store_backed_server['url']}/ga4gh/drs/service-info", timeout=30)
        assert resp.status_code == 200
        data = resp.json()
        assert data["type"]["artifact"] == "drs"


# --- merged from test_versioning_integration.py ------------------------------
# Integration tests for multi-version asset class and recipe support: cross-entity
# resolution, real fixture data, pinned version syntax, overwrite behavior, and
# data channel sync. Basic get/remove/semver-validation logic is covered by unit
# tests in test_recipes.py.


def _write_yaml(path: Path, data: dict) -> Path:
    """Write a YAML file and return its path."""
    path.write_text(yaml.dump(data, default_flow_style=False))
    return path


def _make_asset_class_yaml(name: str, version: str, **kwargs) -> dict:
    """Create an asset class YAML dict."""
    return {
        "name": name,
        "version": version,
        "description": kwargs.get("description", f"{name} v{version}"),
        "serving_modes": kwargs.get("serving_modes", ["file"]),
        "seek_keys": kwargs.get(
            "seek_keys",
            {
                "main": {
                    "value": "{genome}.dat",
                    "description": "Main file",
                    "type": "file",
                }
            },
        ),
    }


def _make_recipe_yaml(
    name: str,
    version: str,
    output_asset_class: str,
    input_assets: dict | None = None,
) -> dict:
    """Create a recipe YAML dict."""
    return {
        "name": name,
        "version": version,
        "description": f"{name} v{version}",
        "output_asset_class": output_asset_class,
        "command_templates": ["echo test"],
        "default_asset": "default",
        "input_params": None,
        "input_files": None,
        "input_assets": input_assets,
        "docker_image": None,
    }


# =============================================================================
# Asset Class Versioning (integration-only scenarios)
# =============================================================================


class TestRecipeWithMultiVersionAssetClass:
    """Recipes that reference multi-version asset classes (integration tier)."""

    def test_real_fasta_multiversion_with_bowtie2_recipe(self, refgenie_instance, tmp_path):
        """Reproduce the original crash: add fasta v0.1.0 + v0.2.0, then add bowtie2 recipe."""
        rg = refgenie_instance

        # Register fasta v0.1.0 (the standard one)
        register_fasta(rg, TESTS_DATA_DIR)

        # Add fasta v0.2.0 as a second version
        fasta_v2 = _write_yaml(
            tmp_path / "fasta_v2.yaml",
            _make_asset_class_yaml(
                "fasta",
                "0.2.0",
                description="FASTA v2",
                serving_modes=["file", "archive"],
                seek_keys={
                    "fasta": {
                        "value": "{genome}.fa",
                        "description": "FASTA file",
                        "type": "file",
                    },
                    "fai": {
                        "value": "{genome}.fa.fai",
                        "description": "FASTA index",
                        "type": "file",
                    },
                    "chrom_sizes": {
                        "value": "{genome}.chrom.sizes",
                        "description": "Chrom sizes",
                        "type": "file",
                    },
                },
            ),
        )
        rg.asset_class.add(fasta_v2)

        # Now add bowtie2 asset class and recipe
        # The bowtie2 recipe references asset_class "fasta" as input
        # This used to crash with MultipleResultsFound
        rg.asset_class.add(TESTS_DATA_DIR / "bowtie2_index_asset_class.yaml")
        rg.recipe.add(TESTS_DATA_DIR / "bowtie2_index_asset_recipe.yaml")

        recipe = rg.recipe.get("bowtie2_index")
        assert recipe is not None
        assert recipe.output_asset_class.name == "bowtie2_index"
