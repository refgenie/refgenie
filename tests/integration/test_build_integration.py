"""
End-to-end integration tests: build -> archive -> serve -> pull -> verify.

Tests the full refgenie workflow using real commands:
- TestRealBuild: Runs actual `refgenie build` subprocess and verifies output
- TestAddArchiveServePull: Full workflow using real build + archive + serve + pull
- TestAddAssetPythonAPI: Uses RefGenConf.add() Python API directly

No backwards compatibility is maintained. This is developmental software.
"""

from __future__ import annotations

import gzip
import hashlib
import os
import shutil
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

# Test genome alias
GENOME_ALIAS = "t7"

# Expected genome digest from building t7.fa.gz (confirmed in CI workflow)
EXPECTED_GENOME_DIGEST = "6c5f19c9c2850e62cc3f89b04047fa05eee911662bd77905"

# Path to the test data directory (relative to repo root)
TESTS_DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
T7_FASTA_GZ = os.path.join(TESTS_DATA_DIR, "t7.fa.gz")
RECIPE_PARENT = os.path.join(TESTS_DATA_DIR, "recipe_parent.json")


def _find_free_port():
    """Find a free TCP port on localhost."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def _read_fasta_content(gz_path):
    """Read and return the decompressed content of a gzipped FASTA file."""
    with gzip.open(gz_path, "rt") as f:
        return f.read()


# =============================================================================
# Fixtures: real build -> archive -> serve
# =============================================================================


@pytest.fixture(scope="session")
def real_build_config(tmp_path_factory):
    """Build a real genome asset using `refgenie build` subprocess.

    Runs `refgenie init` followed by `refgenie build t7/fasta` with the
    test recipe and test FASTA file. This exercises the full build pipeline
    including pypiper, recipe resolution, and config writing.

    Returns a dict with config_path, genome_name, genome_digest, temp_dir,
    and the build return code.
    """
    base = tmp_path_factory.mktemp("real_build")
    config_path = str(base / "genome_config.yaml")

    # refgenie init
    result = subprocess.run(
        [
            sys.executable, "-m", "refgenie", "init",
            "-c", config_path,
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"Init failed: {result.stderr}"

    # refgenie build
    result = subprocess.run(
        [
            sys.executable, "-m", "refgenie", "build",
            "-c", config_path,
            f"{GENOME_ALIAS}/fasta",
            "--files", f"fasta={T7_FASTA_GZ}",
            "--recipe", RECIPE_PARENT,
        ],
        capture_output=True,
        text=True,
    )

    # Get the genome digest dynamically via refgenie id
    id_result = subprocess.run(
        [
            sys.executable, "-m", "refgenie", "id",
            "-c", config_path,
            f"{GENOME_ALIAS}/fasta",
        ],
        capture_output=True,
        text=True,
    )

    # Compute genome digest from refgenie alias get
    alias_result = subprocess.run(
        [
            sys.executable, "-m", "refgenie", "alias", "get",
            "-c", config_path,
            "-a", GENOME_ALIAS,
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
    from yacman import write_lock

    from refgenieserver.server_builder import archive

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


# Separate fixture for add() API tests (not used by build/archive flow)
@pytest.fixture(scope="session")
def added_asset_config(tmp_path_factory):
    """Create a genome asset using RefGenConf.add().

    This exercises the add() code path with write_lock, separate from
    the build/archive flow fixtures. Uses a small inline FASTA for isolation.
    """
    from refgenconf import RefGenConf
    from refgenconf.const import (
        CFG_FOLDER_KEY,
        CFG_GENOMES_KEY,
        CFG_SERVERS_KEY,
        CFG_VERSION_KEY,
    )

    FASTA_CONTENT = """>chrM
GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTT
CGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTC
GCAGTATCTGTCTTTGATTCCTGCCTCATTCTATTATTTATCGCACCTACGTTCAATATT
"""
    FAI_CONTENT = "chrM\t180\t6\t60\t61\n"
    CHROM_SIZES_CONTENT = "chrM\t180\n"

    base = tmp_path_factory.mktemp("add_test")
    genome_folder = base / "genomes"
    genome_folder.mkdir()
    (genome_folder / "data").mkdir()
    (genome_folder / "alias").mkdir()

    input_fasta = base / "input.fa"
    input_fasta.write_text(FASTA_CONTENT)

    config_path = base / "genome_config.yaml"
    config = {
        CFG_VERSION_KEY: 0.4,
        CFG_FOLDER_KEY: str(genome_folder),
        CFG_SERVERS_KEY: [],
        CFG_GENOMES_KEY: None,
    }
    rgc = RefGenConf(entries=config)
    rgc.initialize_config_file(str(config_path))

    rgc = RefGenConf.from_yaml_file(str(config_path))
    genome_digest = rgc.initialize_genome(
        fasta_path=str(input_fasta),
        alias="test_mito",
        fasta_unzipped=True,
    )

    # Create asset directory structure
    asset_dir = genome_folder / "data" / genome_digest / "fasta" / "default"
    asset_dir.mkdir(parents=True, exist_ok=True)

    fa_name = f"{genome_digest}.fa"
    shutil.copy(str(input_fasta), str(asset_dir / fa_name))
    (asset_dir / f"{fa_name}.fai").write_text(FAI_CONTENT)
    (asset_dir / f"{genome_digest}.chrom.sizes").write_text(CHROM_SIZES_CONTENT)

    # Use rgc.add() -- the path must be relative to genome_folder
    rel_path = os.path.join("data", genome_digest, "fasta", "default")
    rgc = RefGenConf.from_yaml_file(str(config_path))
    rgc.add(
        path=rel_path,
        genome="test_mito",
        asset="fasta",
        tag="default",
        seek_keys={
            "fasta": f"{genome_digest}.fa",
            "fai": f"{genome_digest}.fa.fai",
            "chrom_sizes": f"{genome_digest}.chrom.sizes",
        },
        force=True,
    )

    return {
        "config_path": str(config_path),
        "genome_folder": str(genome_folder),
        "genome_digest": genome_digest,
        "genome_alias": "test_mito",
    }


# =============================================================================
# Test group 1: Real build tests
# =============================================================================


class TestRealBuild:
    """Tests that run the actual `refgenie build` pipeline end-to-end.

    These tests use a real `refgenie build` subprocess with the test recipe
    and test FASTA file, then verify the build output is correct.
    """

    def test_real_build_exit_code_zero(self, real_build_config):
        """The build command exits with return code 0."""
        assert real_build_config["build_returncode"] == 0, (
            f"Build failed with return code {real_build_config['build_returncode']}.\n"
            f"stdout: {real_build_config['build_stdout']}\n"
            f"stderr: {real_build_config['build_stderr']}"
        )

    def test_real_build_config_has_genome(self, real_build_config):
        """After build, the config contains the genome digest."""
        from refgenconf import RefGenConf
        from refgenconf.const import CFG_GENOMES_KEY

        rgc = RefGenConf.from_yaml_file(real_build_config["config_path"])
        genome_digest = real_build_config["genome_digest"]
        assert genome_digest in rgc[CFG_GENOMES_KEY], (
            f"Genome digest '{genome_digest}' not found in config genomes: "
            f"{list(rgc[CFG_GENOMES_KEY].keys())}"
        )

    def test_real_build_config_has_fasta_asset(self, real_build_config):
        """After build, the fasta asset exists under the genome in the config."""
        from refgenconf import RefGenConf
        from refgenconf.const import CFG_ASSETS_KEY, CFG_GENOMES_KEY

        rgc = RefGenConf.from_yaml_file(real_build_config["config_path"])
        genome_digest = real_build_config["genome_digest"]
        assert "fasta" in rgc[CFG_GENOMES_KEY][genome_digest][CFG_ASSETS_KEY]

    def test_real_build_seek_returns_valid_path(self, real_build_config):
        """After build, seek returns a path to an existing FASTA file."""
        from refgenconf import RefGenConf

        rgc = RefGenConf.from_yaml_file(real_build_config["config_path"])
        fasta_path = rgc.seek(GENOME_ALIAS, "fasta")
        assert os.path.isfile(fasta_path), f"seek returned nonexistent path: {fasta_path}"

    def test_real_build_fasta_content_correct(self, real_build_config):
        """The built FASTA file contains the expected sequences from t7.fa.gz."""
        from refgenconf import RefGenConf

        rgc = RefGenConf.from_yaml_file(real_build_config["config_path"])
        fasta_path = rgc.seek(GENOME_ALIAS, "fasta")

        with open(fasta_path) as f:
            built_content = f.read()

        # Compare with the original gzipped FASTA content
        original_content = _read_fasta_content(T7_FASTA_GZ)
        assert built_content == original_content

    def test_real_build_asset_digest_present(self, real_build_config):
        """After build, asset_digest is set and non-empty in the config."""
        from refgenconf import RefGenConf
        from refgenconf.const import (
            CFG_ASSETS_KEY,
            CFG_ASSET_CHECKSUM_KEY,
            CFG_ASSET_TAGS_KEY,
            CFG_GENOMES_KEY,
        )

        rgc = RefGenConf.from_yaml_file(real_build_config["config_path"])
        genome_digest = real_build_config["genome_digest"]
        tag_data = rgc[CFG_GENOMES_KEY][genome_digest][CFG_ASSETS_KEY]["fasta"][
            CFG_ASSET_TAGS_KEY
        ]["default"]
        assert CFG_ASSET_CHECKSUM_KEY in tag_data
        assert tag_data[CFG_ASSET_CHECKSUM_KEY] is not None
        assert tag_data[CFG_ASSET_CHECKSUM_KEY] != ""

    def test_real_build_alias_resolves(self, real_build_config):
        """After build, the genome alias resolves to the expected digest."""
        from refgenconf import RefGenConf

        rgc = RefGenConf.from_yaml_file(real_build_config["config_path"])
        digest = rgc.get_genome_alias_digest(GENOME_ALIAS)
        assert digest == EXPECTED_GENOME_DIGEST


# =============================================================================
# Test group 2: Build -> Archive -> Serve -> Pull
# =============================================================================


class TestAddArchiveServePull:
    """True end-to-end tests: build -> archive -> serve -> pull -> verify.

    Uses real `refgenie build` output, refgenieserver archive() to create
    archives, create_app() to serve, and RefGenConf.pull() to pull.
    Verifies content matches original.
    """

    def test_server_starts_with_archived_assets(self, build_test_server):
        """The server starts and serves the archived genome."""
        import requests

        resp = requests.get(f"{build_test_server['url']}/v3/genomes/list")
        assert resp.status_code == 200
        data = resp.json()
        assert build_test_server["genome_digest"] in data

    def test_server_serves_asset_archive(self, build_test_server):
        """The server serves the archived asset for download."""
        import requests

        resp = requests.get(
            f"{build_test_server['url']}/v3/assets/archive/"
            f"{build_test_server['genome_digest']}/fasta",
            params={"tag": "default"},
        )
        assert resp.status_code == 200
        # Verify it's a valid gzip file
        assert resp.content[:2] == b"\x1f\x8b"

    def test_pull_fasta_content_matches_original(self, build_test_server, tmp_path):
        """After pull, FASTA content matches the originally built file."""
        from refgenconf import RefGenConf
        from refgenconf.const import CFG_FOLDER_KEY, CFG_SERVERS_KEY, CFG_VERSION_KEY

        genome_folder = tmp_path / "pulled_genomes"
        genome_folder.mkdir()
        (genome_folder / "data").mkdir()
        (genome_folder / "alias").mkdir()

        config_path = tmp_path / "client.yaml"
        config = {
            CFG_VERSION_KEY: 0.4,
            CFG_FOLDER_KEY: str(genome_folder),
            CFG_SERVERS_KEY: [build_test_server["url"]],
            "genomes": {},
        }
        with open(str(config_path), "w") as f:
            yaml.dump(config, f, default_flow_style=False)

        rgc = RefGenConf.from_yaml_file(str(config_path))
        result = rgc.pull(GENOME_ALIAS, "fasta", "default", force=True)
        assert result is not None

        # Verify FASTA content matches original
        fasta_path = rgc.seek(GENOME_ALIAS, "fasta")
        with open(fasta_path) as f:
            pulled_content = f.read()

        original_content = _read_fasta_content(T7_FASTA_GZ)
        assert pulled_content == original_content

    def test_archive_checksum_matches_download(self, build_test_server):
        """The archive checksum reported by server matches actual download."""
        import requests

        genome_digest = build_test_server["genome_digest"]

        # Get reported digest
        resp_digest = requests.get(
            f"{build_test_server['url']}/v3/assets/archive_digest/"
            f"{genome_digest}/fasta"
        )
        reported_digest = resp_digest.text

        # Download archive and compute digest
        resp_archive = requests.get(
            f"{build_test_server['url']}/v3/assets/archive/"
            f"{genome_digest}/fasta",
            params={"tag": "default"},
        )
        computed_digest = hashlib.md5(resp_archive.content).hexdigest()
        assert computed_digest == reported_digest


# =============================================================================
# Test group 3: Full workflow via CLI commands
# =============================================================================


class TestFullWorkflowCLI:
    """Test the complete workflow via CLI: init -> pull -> seek -> verify."""

    def test_cli_pull_and_seek(self, build_test_server, tmp_path):
        """refgenie pull + seek via CLI returns correct file content."""
        client_cfg = str(tmp_path / "client.yaml")
        genome_folder = str(tmp_path / "pulled")

        # Init
        result = subprocess.run(
            [
                sys.executable, "-m", "refgenie", "init",
                "-c", client_cfg, "-g", genome_folder,
                "-s", build_test_server["url"],
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, f"Init failed: {result.stderr}"

        # Pull
        result = subprocess.run(
            [
                sys.executable, "-m", "refgenie", "pull",
                "-c", client_cfg,
                f"{GENOME_ALIAS}/fasta:default",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, f"Pull failed: {result.stderr}"

        # Seek
        result = subprocess.run(
            [
                sys.executable, "-m", "refgenie", "seek",
                "-c", client_cfg,
                f"{GENOME_ALIAS}/fasta",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, f"Seek failed: {result.stderr}"
        fasta_path = result.stdout.strip()
        assert os.path.isfile(fasta_path), (
            f"Seek returned nonexistent path: {fasta_path}"
        )

        # Verify content
        with open(fasta_path) as f:
            content = f.read()
        original_content = _read_fasta_content(T7_FASTA_GZ)
        assert content == original_content

    def test_cli_listr_shows_built_asset(self, build_test_server, tmp_path):
        """refgenie listr shows the asset that was built and archived."""
        client_cfg = str(tmp_path / "client.yaml")
        genome_folder = str(tmp_path / "listed")

        subprocess.run(
            [
                sys.executable, "-m", "refgenie", "init",
                "-c", client_cfg, "-g", genome_folder,
                "-s", build_test_server["url"],
            ],
            capture_output=True, text=True, check=True,
        )

        result = subprocess.run(
            [
                sys.executable, "-m", "refgenie", "listr",
                "-c", client_cfg,
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, f"listr failed: {result.stderr}"
        assert (
            GENOME_ALIAS in result.stdout
            or build_test_server["genome_digest"] in result.stdout
        )

    def test_cli_pull_and_list_shows_local_asset(self, build_test_server, tmp_path):
        """After CLI pull, refgenie list shows the asset locally."""
        client_cfg = str(tmp_path / "client.yaml")
        genome_folder = str(tmp_path / "local")

        subprocess.run(
            [
                sys.executable, "-m", "refgenie", "init",
                "-c", client_cfg, "-g", genome_folder,
                "-s", build_test_server["url"],
            ],
            capture_output=True, text=True, check=True,
        )

        subprocess.run(
            [
                sys.executable, "-m", "refgenie", "pull",
                "-c", client_cfg,
                f"{GENOME_ALIAS}/fasta:default",
            ],
            capture_output=True, text=True, check=True,
        )

        result = subprocess.run(
            [
                sys.executable, "-m", "refgenie", "list",
                "-c", client_cfg,
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "fasta" in result.stdout


# =============================================================================
# Test group 4: Asset registration via add() Python API
# =============================================================================


class TestAddAssetPythonAPI:
    """Test that RefGenConf.add() correctly registers assets.

    These tests verify the add() code path which uses write_lock.
    Uses a small inline FASTA for isolation from the build tests.
    """

    def test_add_creates_config_entry(self, added_asset_config):
        """After add(), the config contains the asset entry."""
        from refgenconf import RefGenConf
        from refgenconf.const import CFG_ASSETS_KEY, CFG_GENOMES_KEY

        rgc = RefGenConf.from_yaml_file(added_asset_config["config_path"])
        genome_digest = added_asset_config["genome_digest"]

        assert genome_digest in rgc[CFG_GENOMES_KEY]
        assert "fasta" in rgc[CFG_GENOMES_KEY][genome_digest][CFG_ASSETS_KEY]

    def test_add_registers_seek_keys(self, added_asset_config):
        """After add(), seek keys are registered correctly."""
        from refgenconf import RefGenConf

        rgc = RefGenConf.from_yaml_file(added_asset_config["config_path"])
        alias = added_asset_config["genome_alias"]

        fasta_path = rgc.seek(alias, "fasta", seek_key="fasta")
        assert os.path.exists(fasta_path)

        fai_path = rgc.seek(alias, "fasta", seek_key="fai")
        assert os.path.exists(fai_path)

        sizes_path = rgc.seek(alias, "fasta", seek_key="chrom_sizes")
        assert os.path.exists(sizes_path)

    def test_add_preserves_genome_alias(self, added_asset_config):
        """After add(), the genome alias resolves to the correct digest."""
        from refgenconf import RefGenConf

        rgc = RefGenConf.from_yaml_file(added_asset_config["config_path"])
        alias = added_asset_config["genome_alias"]
        digest = rgc.get_genome_alias_digest(alias=alias)
        assert digest == added_asset_config["genome_digest"]

    def test_add_computes_asset_digest(self, added_asset_config):
        """After add(), the asset has a computed digest."""
        from refgenconf import RefGenConf
        from refgenconf.const import (
            CFG_ASSETS_KEY,
            CFG_ASSET_CHECKSUM_KEY,
            CFG_ASSET_TAGS_KEY,
            CFG_GENOMES_KEY,
        )

        rgc = RefGenConf.from_yaml_file(added_asset_config["config_path"])
        genome_digest = added_asset_config["genome_digest"]

        tag_data = rgc[CFG_GENOMES_KEY][genome_digest][CFG_ASSETS_KEY]["fasta"][
            CFG_ASSET_TAGS_KEY
        ]["default"]
        assert CFG_ASSET_CHECKSUM_KEY in tag_data
        assert tag_data[CFG_ASSET_CHECKSUM_KEY] != ""
        assert tag_data[CFG_ASSET_CHECKSUM_KEY] is not None


# =============================================================================
# Test group 5: Archive integrity tests
# =============================================================================


class TestArchiveCreation:
    """Verify that refgenieserver archive produces correct archives."""

    def test_archive_file_exists(self, archived_asset_config):
        """The archive .tgz file was created by archive()."""
        genome_digest = archived_asset_config["genome_digest"]
        tgz_path = os.path.join(
            archived_asset_config["archive_dir"], genome_digest, "fasta__default.tgz"
        )
        assert os.path.exists(tgz_path)
        assert os.path.getsize(tgz_path) > 0

    def test_archive_is_valid_tarball(self, archived_asset_config):
        """The archive is a valid gzipped tarball."""
        import tarfile

        genome_digest = archived_asset_config["genome_digest"]
        tgz_path = os.path.join(
            archived_asset_config["archive_dir"], genome_digest, "fasta__default.tgz"
        )
        with tarfile.open(tgz_path, "r:gz") as tar:
            names = tar.getnames()
        assert len(names) > 0
        assert any("default" in n for n in names)

    def test_archive_contains_fasta(self, archived_asset_config):
        """The archive contains the FASTA file."""
        import tarfile

        genome_digest = archived_asset_config["genome_digest"]
        tgz_path = os.path.join(
            archived_asset_config["archive_dir"], genome_digest, "fasta__default.tgz"
        )
        with tarfile.open(tgz_path, "r:gz") as tar:
            names = tar.getnames()
        assert any(name.endswith(".fa") for name in names)

    def test_server_config_has_archive_metadata(self, archived_asset_config):
        """The server config produced by archive() has checksum and size."""
        from refgenconf.const import (
            CFG_ARCHIVE_CHECKSUM_KEY,
            CFG_ARCHIVE_SIZE_KEY,
            CFG_ASSETS_KEY,
            CFG_ASSET_TAGS_KEY,
            CFG_GENOMES_KEY,
        )

        with open(archived_asset_config["server_config_path"]) as f:
            config = yaml.safe_load(f)

        genome_digest = archived_asset_config["genome_digest"]
        tag_data = config[CFG_GENOMES_KEY][genome_digest][CFG_ASSETS_KEY]["fasta"][
            CFG_ASSET_TAGS_KEY
        ]["default"]

        assert CFG_ARCHIVE_CHECKSUM_KEY in tag_data
        assert tag_data[CFG_ARCHIVE_CHECKSUM_KEY] is not None
        assert len(tag_data[CFG_ARCHIVE_CHECKSUM_KEY]) == 32  # md5 hex

        assert CFG_ARCHIVE_SIZE_KEY in tag_data
        assert tag_data[CFG_ARCHIVE_SIZE_KEY] is not None
