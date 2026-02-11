"""
Cross-package integration test fixtures for refgenie + refgenconf + refgenieserver.

Creates a minimal FastAPI server mimicking refgenieserver's API, serves archived
genome assets, and provides a RefGenConf client subscribed to the test server.

Uses a real threaded uvicorn server so the refgenconf client exercises its actual
HTTP code paths (requests.get, urlopen, urlretrieve).

Prerequisites:
    pip install -e ../refgenconf
    pip install -e ../refgenieserver
    pip install -e .

Or use the runner script:
    ./tests/scripts/test-integration.sh
"""

from __future__ import annotations

import hashlib
import logging
import os
import socket
import tarfile
import threading
import time
from copy import deepcopy
from pathlib import Path
from typing import Any

import pytest
import yaml

# Skip all integration tests unless explicitly enabled
pytestmark = pytest.mark.skipif(
    os.getenv("RUN_INTEGRATION_TESTS") != "true",
    reason="Integration tests disabled. Run ./tests/scripts/test-integration.sh",
)

# Test genome: rCRSd (revised Cambridge Reference Sequence, doubled)
# This digest is computed by SeqColClient.load_fasta() from the FASTA_CONTENT below.
# It MUST match what initialize_genome() computes, or pull tests will fail.
GENOME_DIGEST = "1d5d914dc1bfbeb42f568f4fd6fe7af14d30f38860939f69"
GENOME_ALIAS = "rCRSd"
GENOME_DESC = "The revised cambridge reference sequence"

# Minimal FASTA content for testing
FASTA_CONTENT = """>chrM
GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTT
CGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTC
GCAGTATCTGTCTTTGATTCCTGCCTCATTCTATTATTTATCGCACCTACGTTCAATATT
"""


@pytest.fixture(scope="session")
def archived_genome_dir(tmp_path_factory):
    """Create a temp directory with an archived genome asset (.tgz file).

    Simulates what `refgenieserver archive` produces: a .tgz archive of the
    asset directory with computed checksums.

    Returns a dict with 'path', 'archive_checksum', and 'archive_size'.
    """
    base = tmp_path_factory.mktemp("server_genomes")

    # Create the asset directory structure:
    # {genome_digest}/fasta/default/
    asset_dir = base / GENOME_DIGEST / "fasta" / "default"
    asset_dir.mkdir(parents=True)

    # Create minimal FASTA files
    fa_name = f"{GENOME_DIGEST}.fa"
    fai_name = f"{GENOME_DIGEST}.fa.fai"
    chrom_sizes_name = f"{GENOME_DIGEST}.chrom.sizes"

    fa_path = asset_dir / fa_name
    fa_path.write_text(FASTA_CONTENT)

    # Create a minimal .fai index
    fai_path = asset_dir / fai_name
    fai_path.write_text("chrM\t180\t6\t60\t61\n")

    # Create chrom sizes
    chrom_sizes_path = asset_dir / chrom_sizes_name
    chrom_sizes_path.write_text("chrM\t180\n")

    # Create the .tgz archive (what refgenieserver serves for download)
    archive_dir = base / GENOME_DIGEST
    tgz_path = archive_dir / "fasta__default.tgz"

    # The tar should contain the "default" directory
    with tarfile.open(str(tgz_path), "w:gz") as tar:
        tar.add(str(asset_dir), arcname="default")

    # Compute archive checksum (md5) — matches how refgenieserver does it
    md5 = hashlib.md5()
    with open(str(tgz_path), "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            md5.update(chunk)
    archive_checksum = md5.hexdigest()
    archive_size = os.path.getsize(str(tgz_path))

    return {
        "path": base,
        "archive_checksum": archive_checksum,
        "archive_size": str(archive_size),
    }


@pytest.fixture(scope="session")
def test_server_app(archived_genome_dir):
    """Create a minimal FastAPI app mimicking refgenieserver's key API endpoints.

    Uses the same operation IDs as refgenieserver so that refgenconf's
    `construct_request_url()` can discover endpoints via /openapi.json.
    """
    from fastapi import FastAPI, Query
    from fastapi.responses import FileResponse, Response

    from refgenconf.const import (
        API_ID_ALIASES_DICT,
        API_ID_ARCHIVE,
        API_ID_ASSET_ATTRS,
        API_ID_ASSETS,
        API_ID_DEFAULT_TAG,
        API_ID_DIGEST,
        API_ID_GENOMES_DICT,
        API_ID_GENOME_ATTRS,
        API_VERSION,
        PRIVATE_API,
    )

    app = FastAPI(title="Test refgenieserver", version="0.0.1")

    genome_dir = archived_genome_dir["path"]
    archive_checksum = archived_genome_dir["archive_checksum"]
    archive_size = archived_genome_dir["archive_size"]

    # Test data matching genomes.yaml structure
    asset_data = {
        "seek_keys": {
            "fasta": f"{GENOME_DIGEST}.fa",
            "fai": f"{GENOME_DIGEST}.fa.fai",
            "chrom_sizes": f"{GENOME_DIGEST}.chrom.sizes",
            "dir": ".",
        },
        "asset_parents": [],
        "asset_path": "fasta",
        "asset_digest": "4eb430296bc02ed7e4006624f1d5ac53",
        "archive_digest": archive_checksum,
        "archive_size": archive_size,
    }

    genome_attrs = {
        "genome_description": GENOME_DESC,
        "aliases": [GENOME_ALIAS],
    }

    # --- API endpoints matching refgenieserver's operation IDs ---

    @app.get(
        "/v3/genomes/list",
        operation_id=API_VERSION + "custom_Id_genomes_list",
    )
    async def list_genomes():
        return [GENOME_DIGEST]

    @app.get(
        "/v3/genomes/alias_dict",
        operation_id=API_VERSION + API_ID_ALIASES_DICT,
    )
    async def get_alias_dict():
        return {GENOME_DIGEST: [GENOME_ALIAS]}

    @app.get(
        "/v3/assets/list",
        operation_id=API_VERSION + API_ID_ASSETS,
    )
    async def list_assets(includeSeekKeys: bool = Query(False)):
        if includeSeekKeys:
            return {GENOME_DIGEST: {"fasta": {"default": asset_data["seek_keys"]}}}
        return {GENOME_DIGEST: ["fasta"]}

    @app.get(
        "/v3/assets/attrs/{genome}/{asset}",
        operation_id=API_VERSION + API_ID_ASSET_ATTRS,
    )
    async def get_asset_attrs(genome: str, asset: str, tag: str = Query(None)):
        if genome == GENOME_DIGEST and asset == "fasta":
            return asset_data
        return Response(status_code=404)

    @app.get(
        "/v3/genomes/attrs/{genome}",
        operation_id=API_VERSION + API_ID_GENOME_ATTRS,
    )
    async def get_genome_attrs(genome: str):
        if genome == GENOME_DIGEST:
            return genome_attrs
        return Response(status_code=404)

    @app.get(
        "/v3/assets/default_tag/{genome}/{asset}",
        operation_id=API_VERSION + API_ID_DEFAULT_TAG,
    )
    async def get_default_tag(genome: str, asset: str):
        return Response(content="default", media_type="text/plain")

    @app.get(
        "/v3/assets/asset_digest/{genome}/{asset}",
        operation_id=API_VERSION + API_ID_DIGEST,
    )
    async def get_asset_digest(genome: str, asset: str, tag: str = Query(None)):
        if genome == GENOME_DIGEST and asset == "fasta":
            return Response(
                content=asset_data["asset_digest"], media_type="text/plain"
            )
        return Response(status_code=404)

    @app.get(
        "/v3/assets/archive/{genome}/{asset}",
        operation_id=API_VERSION + API_ID_ARCHIVE,
    )
    async def download_archive(genome: str, asset: str, tag: str = Query(None)):
        if genome == GENOME_DIGEST and asset == "fasta":
            tgz_path = str(genome_dir / GENOME_DIGEST / "fasta__default.tgz")
            return FileResponse(
                tgz_path,
                filename="fasta__default.tgz",
                media_type="application/octet-stream",
            )
        return Response(status_code=404)

    @app.get(
        "/v3/genomes/genome_digest/{alias}",
        operation_id=API_VERSION + "custom_Id_alias_digest",
    )
    async def get_genome_digest(alias: str):
        if alias == GENOME_ALIAS or alias == GENOME_DIGEST:
            return Response(content=GENOME_DIGEST, media_type="text/plain")
        return Response(status_code=404)

    # --- Private API endpoint (used by CLI listr command) ---

    @app.get(
        "/genomes/dict",
        operation_id=PRIVATE_API + API_ID_GENOMES_DICT,
    )
    async def get_genomes_dict():
        return {
            GENOME_DIGEST: {
                "aliases": [GENOME_ALIAS],
                "genome_description": GENOME_DESC,
                "assets": {
                    "fasta": {
                        "default_tag": "default",
                        "tags": {
                            "default": {
                                "asset_path": "fasta",
                                "asset_digest": "4eb430296bc02ed7e4006624f1d5ac53",
                                "archive_digest": archive_checksum,
                                "archive_size": archive_size,
                                "seek_keys": asset_data["seek_keys"],
                                "asset_parents": [],
                            }
                        },
                    }
                },
            }
        }

    return app


def _find_free_port():
    """Find a free TCP port on localhost."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


@pytest.fixture(scope="session")
def test_server_url(test_server_app):
    """Start the test server on a random port in a background thread.

    Yields the server URL (e.g. http://127.0.0.1:12345).
    """
    import uvicorn

    port = _find_free_port()
    url = f"http://127.0.0.1:{port}"

    config = uvicorn.Config(
        test_server_app,
        host="127.0.0.1",
        port=port,
        log_level="warning",
    )
    server = uvicorn.Server(config)

    thread = threading.Thread(target=server.run, daemon=True)
    thread.start()

    # Wait for server to be ready
    for _ in range(50):
        try:
            with socket.create_connection(("127.0.0.1", port), timeout=0.1):
                break
        except (ConnectionRefusedError, OSError):
            time.sleep(0.1)
    else:
        raise RuntimeError(f"Test server failed to start on port {port}")

    yield url

    server.should_exit = True
    thread.join(timeout=5)


@pytest.fixture(scope="session")
def server_config_path(archived_genome_dir, tmp_path_factory):
    """Create a server config YAML file (what refgenieserver would use)."""
    config_dir = tmp_path_factory.mktemp("server_config")
    config_path = config_dir / "server_config.yaml"

    archive_checksum = archived_genome_dir["archive_checksum"]
    archive_size = archived_genome_dir["archive_size"]

    config = {
        "config_version": 0.4,
        "genome_folder": str(archived_genome_dir["path"]),
        "genome_servers": [],
        "genomes": {
            GENOME_DIGEST: {
                "aliases": [GENOME_ALIAS],
                "genome_description": GENOME_DESC,
                "assets": {
                    "fasta": {
                        "default_tag": "default",
                        "tags": {
                            "default": {
                                "asset_path": "fasta",
                                "asset_digest": "4eb430296bc02ed7e4006624f1d5ac53",
                                "archive_digest": archive_checksum,
                                "archive_size": archive_size,
                                "seek_keys": {
                                    "fasta": f"{GENOME_DIGEST}.fa",
                                    "fai": f"{GENOME_DIGEST}.fa.fai",
                                    "chrom_sizes": f"{GENOME_DIGEST}.chrom.sizes",
                                    "dir": ".",
                                },
                                "asset_parents": [],
                            }
                        },
                    }
                },
            }
        },
    }

    with open(str(config_path), "w") as f:
        yaml.dump(config, f, default_flow_style=False)

    return str(config_path)


@pytest.fixture
def client_rgc(test_server_url, tmp_path):
    """Create a writable RefGenConf client subscribed to the test server.

    This is the main fixture for testing client-server operations.
    """
    from refgenconf import RefGenConf
    from refgenconf.const import CFG_SERVERS_KEY, CFG_VERSION_KEY, CFG_FOLDER_KEY

    genome_folder = tmp_path / "client_genomes"
    genome_folder.mkdir()
    (genome_folder / "data").mkdir()
    (genome_folder / "alias").mkdir()

    config_path = tmp_path / "refgenie.yaml"
    config = {
        CFG_VERSION_KEY: 0.4,
        CFG_FOLDER_KEY: str(genome_folder),
        CFG_SERVERS_KEY: [test_server_url],
        "genomes": {},
    }

    with open(str(config_path), "w") as f:
        yaml.dump(config, f, default_flow_style=False)

    return RefGenConf.from_yaml_file(str(config_path))


@pytest.fixture
def client_config_path(test_server_url, tmp_path):
    """Create a client config file and return its path (for CLI tests)."""
    from refgenconf.const import CFG_SERVERS_KEY, CFG_VERSION_KEY, CFG_FOLDER_KEY

    genome_folder = tmp_path / "cli_genomes"
    genome_folder.mkdir()
    (genome_folder / "data").mkdir()
    (genome_folder / "alias").mkdir()

    config_path = tmp_path / "refgenie.yaml"
    config = {
        CFG_VERSION_KEY: 0.4,
        CFG_FOLDER_KEY: str(genome_folder),
        CFG_SERVERS_KEY: [test_server_url],
        "genomes": {},
    }

    with open(str(config_path), "w") as f:
        yaml.dump(config, f, default_flow_style=False)

    return str(config_path)
