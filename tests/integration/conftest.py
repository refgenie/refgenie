"""
Cross-package integration test fixtures for refgenie + refgenconf + refgenieserver.

Uses real `refgenie build` and real `archive()` instead of hand-crafted test data.
This ensures integration tests catch real bugs in the build pipeline, archive
creation, server startup, and client interaction.

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

import gzip
import os
import socket
import subprocess
import sys
import threading
import time

import pytest
import yaml

# Skip all integration tests unless explicitly enabled
pytestmark = pytest.mark.skipif(
    os.getenv("RUN_INTEGRATION_TESTS") != "true",
    reason="Integration tests disabled. Run ./tests/scripts/test-integration.sh",
)

# Test genome alias and description
GENOME_ALIAS = "rCRSd"
GENOME_DESC = "The revised cambridge reference sequence"

# Path to the test data directory
TESTS_DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
T7_FASTA_GZ = os.path.join(TESTS_DATA_DIR, "t7.fa.gz")
RECIPE_PARENT = os.path.join(TESTS_DATA_DIR, "recipe_parent.json")


def _read_fasta_content(gz_path):
    """Read and return the decompressed content of a gzipped FASTA file."""
    with gzip.open(gz_path, "rt") as f:
        return f.read()


def _find_free_port():
    """Find a free TCP port on localhost."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


# GENOME_DIGEST is deterministic for a given FASTA file. This digest is
# computed by initialize_genome() from the t7.fa.gz test data and is also
# confirmed by the CI workflow in test-refgenie-cli.yml.
GENOME_DIGEST = "6c5f19c9c2850e62cc3f89b04047fa05eee911662bd77905"

# FASTA content read from the test file at import time
FASTA_CONTENT = _read_fasta_content(T7_FASTA_GZ)


# =============================================================================
# Fixtures: real build -> archive -> serve (shared by all test files)
# =============================================================================


@pytest.fixture(scope="session")
def real_build_config(tmp_path_factory):
    """Build a real genome asset using `refgenie build` subprocess.

    Runs `refgenie init` followed by `refgenie build t7/fasta` with the
    test recipe and test FASTA file. This exercises the full build pipeline
    including pypiper, recipe resolution, and config writing.

    Returns a dict with config_path, genome_name, genome_digest, temp_dir,
    and the build result.
    """
    base = tmp_path_factory.mktemp("real_build")
    config_path = str(base / "genome_config.yaml")

    # refgenie init
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "refgenie",
            "init",
            "-c",
            config_path,
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"Init failed: {result.stderr}"

    # refgenie build
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "refgenie",
            "build",
            "-c",
            config_path,
            f"{GENOME_ALIAS}/fasta",
            "--files",
            f"fasta={T7_FASTA_GZ}",
            "--recipe",
            RECIPE_PARENT,
            "--genome-description",
            GENOME_DESC,
        ],
        capture_output=True,
        text=True,
    )

    # Get the genome digest dynamically via refgenie id
    id_result = subprocess.run(
        [
            sys.executable,
            "-m",
            "refgenie",
            "id",
            "-c",
            config_path,
            f"{GENOME_ALIAS}/fasta",
        ],
        capture_output=True,
        text=True,
    )

    # Compute genome digest from refgenie alias get
    alias_result = subprocess.run(
        [
            sys.executable,
            "-m",
            "refgenie",
            "alias",
            "get",
            "-c",
            config_path,
            "-a",
            GENOME_ALIAS,
        ],
        capture_output=True,
        text=True,
    )
    genome_digest = alias_result.stdout.strip() if alias_result.returncode == 0 else ""

    return {
        "config_path": config_path,
        "genome_name": GENOME_ALIAS,
        "genome_digest": genome_digest,
        "temp_dir": str(base),
        "build_returncode": result.returncode,
        "build_stdout": result.stdout,
        "build_stderr": result.stderr,
        "asset_digest": id_result.stdout.strip() if id_result.returncode == 0 else "",
    }


@pytest.fixture(scope="session")
def archived_asset_config(real_build_config, tmp_path_factory):
    """Archive the built asset using refgenieserver archive.

    Calls the archive() function from refgenieserver.server_builder directly,
    which exercises the real archiving code path including tar creation,
    checksum computation, and config updates.

    Returns dict with server_config_path, archive_dir, and genome_digest.
    """
    from refgenconf import RefGenConf
    from refgenconf.const import CFG_ARCHIVE_KEY
    from refgenieserver.server_builder import archive
    from yacman import write_lock

    assert real_build_config["build_returncode"] == 0, (
        f"Build failed, cannot archive: {real_build_config['build_stderr']}"
    )

    archive_dir = tmp_path_factory.mktemp("archives")

    # Update the config to include archive folder
    rgc = RefGenConf.from_yaml_file(real_build_config["config_path"])
    with write_lock(rgc) as r:
        r[CFG_ARCHIVE_KEY] = str(archive_dir)
        r.write()

    # Run refgenieserver archive
    rgc = RefGenConf.from_yaml_file(real_build_config["config_path"])
    archive(
        rgc=rgc,
        registry_paths=None,  # archive all
        force=True,
        remove=False,
        cfg_path=real_build_config["config_path"],
        genomes_desc=None,
    )

    # Verify the archive was created
    genome_digest = real_build_config["genome_digest"]
    tgz_path = os.path.join(str(archive_dir), genome_digest, "fasta__default.tgz")
    assert os.path.exists(tgz_path), f"Archive not created at {tgz_path}"

    # The server config is written by archive() into the archive dir
    server_config_path = os.path.join(
        str(archive_dir), os.path.basename(real_build_config["config_path"])
    )
    assert os.path.exists(server_config_path), (
        f"Server config not found at {server_config_path}"
    )

    return {
        "server_config_path": server_config_path,
        "archive_dir": str(archive_dir),
        "genome_digest": genome_digest,
    }


@pytest.fixture(scope="session")
def build_test_server(archived_asset_config):
    """Start a real refgenieserver serving the archived assets.

    Uses create_app() from refgenieserver.app_factory to create
    a real FastAPI app, then starts it on a random port via uvicorn
    in a background thread.

    Yields dict with url and genome_digest.
    """
    import uvicorn
    from refgenieserver.app_factory import create_app

    app = create_app(
        archived_asset_config["server_config_path"],
        archive_base_dir=archived_asset_config["archive_dir"],
    )

    port = _find_free_port()
    server_url = f"http://127.0.0.1:{port}"

    config = uvicorn.Config(
        app,
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
        raise RuntimeError(f"Build test server failed to start on port {port}")

    yield {
        "url": server_url,
        "genome_digest": archived_asset_config["genome_digest"],
    }

    server.should_exit = True
    thread.join(timeout=5)


# =============================================================================
# Fixtures for test_server_client.py (thin wrappers around shared fixtures)
# =============================================================================


@pytest.fixture(scope="session")
def test_server_url(build_test_server):
    """Return the URL of the shared test server."""
    return build_test_server["url"]


@pytest.fixture(scope="session")
def server_rgc(archived_asset_config):
    """Return the RefGenConf object that the server uses (for testing helpers)."""
    from refgenconf import RefGenConf

    return RefGenConf.from_yaml_file(archived_asset_config["server_config_path"])


@pytest.fixture
def client_rgc(test_server_url, tmp_path):
    """Create a writable RefGenConf client subscribed to the test server.

    This is the main fixture for testing client-server operations.
    """
    from refgenconf import RefGenConf
    from refgenconf.const import CFG_FOLDER_KEY, CFG_SERVERS_KEY, CFG_VERSION_KEY

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
    from refgenconf.const import CFG_FOLDER_KEY, CFG_SERVERS_KEY, CFG_VERSION_KEY

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
