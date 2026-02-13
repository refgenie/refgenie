"""
Integration tests for refgenie build --map / --reduce mode.

Tests the map/reduce build workflow:
- Map step: builds assets and stores metadata in separate per-asset configs
- Reduce step: merges per-asset configs into the master config

No backwards compatibility is maintained. This is developmental software.
"""

from __future__ import annotations

import os
import subprocess
import sys

import pytest

# Skip all integration tests unless explicitly enabled
pytestmark = pytest.mark.skipif(
    os.getenv("RUN_INTEGRATION_TESTS") != "true",
    reason="Integration tests disabled. Run ./tests/scripts/test-integration.sh",
)

# Test genome alias
GENOME_ALIAS = "t7"

# Expected genome digest from building t7.fa.gz
EXPECTED_GENOME_DIGEST = "6c5f19c9c2850e62cc3f89b04047fa05eee911662bd77905"

# Path to the test data directory
TESTS_DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
T7_FASTA_GZ = os.path.join(TESTS_DATA_DIR, "t7.fa.gz")
RECIPE_PARENT = os.path.join(TESTS_DATA_DIR, "recipe_parent.json")


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture(scope="session")
def map_build_config(tmp_path_factory):
    """Build a genome asset using --map mode.

    Runs `refgenie init` followed by `refgenie build --map` with the test
    recipe and test FASTA file. Uses --preserve-map-configs to keep the
    map config file around for inspection.

    Returns a dict with config_path, genome_name, temp_dir, and build results.
    """
    base = tmp_path_factory.mktemp("map_build")
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

    # refgenie build --map --preserve-map-configs
    result = subprocess.run(
        [
            sys.executable, "-m", "refgenie", "build",
            "-c", config_path,
            f"{GENOME_ALIAS}/fasta",
            "--map",
            "--preserve-map-configs",
            "--files", f"fasta={T7_FASTA_GZ}",
            "--recipe", RECIPE_PARENT,
        ],
        capture_output=True,
        text=True,
    )

    return {
        "config_path": config_path,
        "genome_name": GENOME_ALIAS,
        "temp_dir": str(base),
        "build_returncode": result.returncode,
        "build_stdout": result.stdout,
        "build_stderr": result.stderr,
    }


@pytest.fixture(scope="session")
def reduced_config(map_build_config, tmp_path_factory):
    """Run the reduce step on the map build output.

    Creates a copy of the map config and runs --reduce on it.
    Returns a dict with config_path and reduce results.
    """
    import shutil

    # Copy the config so we can run reduce independently
    base = tmp_path_factory.mktemp("reduce_build")
    src_config = map_build_config["config_path"]

    # Create the reduce config by copying the original
    reduce_config_path = str(base / "genome_config.yaml")
    shutil.copy2(src_config, reduce_config_path)

    # Update the genome_folder in the reduce config to point to the real data
    # (the data is in the map build temp dir)
    import yaml

    with open(reduce_config_path) as f:
        config = yaml.safe_load(f)
    # Keep the original genome_folder path since the data is there
    with open(reduce_config_path, "w") as f:
        yaml.dump(config, f, default_flow_style=False)

    # Run reduce
    result = subprocess.run(
        [
            sys.executable, "-m", "refgenie", "build",
            "--reduce",
            "--preserve-map-configs",
            "-c", reduce_config_path,
        ],
        capture_output=True,
        text=True,
    )

    return {
        "config_path": reduce_config_path,
        "reduce_returncode": result.returncode,
        "reduce_stdout": result.stdout,
        "reduce_stderr": result.stderr,
        "original_map_config_path": src_config,
    }


# =============================================================================
# Test group 1: Map step tests
# =============================================================================


class TestMapBuild:
    """Tests for the --map build step."""

    def test_map_build_exit_code_zero(self, map_build_config):
        """The --map build command exits with return code 0."""
        assert map_build_config["build_returncode"] == 0, (
            f"Map build failed with return code {map_build_config['build_returncode']}.\n"
            f"stdout: {map_build_config['build_stdout']}\n"
            f"stderr: {map_build_config['build_stderr']}"
        )

    def test_map_build_creates_map_config(self, map_build_config):
        """After --map build, a map config file exists in the build stats directory."""
        from refgenconf import RefGenConf
        from refgenconf.const import BUILD_MAP_CFG

        rgc = RefGenConf.from_yaml_file(map_build_config["config_path"])
        data_dir = rgc.data_dir

        # The map config should be at:
        # <data_dir>/<genome_digest>/fasta/default/_refgenie_build/_map_build.yaml
        # We need to find it by searching for the genome digest directory
        found_map_config = False
        if os.path.exists(data_dir):
            for genome_dir in os.listdir(data_dir):
                map_cfg_path = os.path.join(
                    data_dir, genome_dir, "fasta", "default",
                    "_refgenie_build", BUILD_MAP_CFG,
                )
                if os.path.exists(map_cfg_path):
                    found_map_config = True
                    break

        assert found_map_config, (
            f"Map config file ({BUILD_MAP_CFG}) not found under {data_dir}"
        )

    def test_map_build_does_not_update_master_config(self, map_build_config):
        """After --map build, the master config has no genome entries."""
        from refgenconf import RefGenConf
        from refgenconf.const import CFG_GENOMES_KEY

        rgc = RefGenConf.from_yaml_file(map_build_config["config_path"])
        # The map step writes to a separate config, not the master
        genomes = rgc.get(CFG_GENOMES_KEY)
        assert genomes is None or len(genomes) == 0, (
            f"Master config should have no genomes after map step, "
            f"but found: {genomes}"
        )

    def test_map_build_map_config_has_genome(self, map_build_config):
        """The map config contains the genome digest and fasta asset."""
        from refgenconf import RefGenConf
        from refgenconf.const import BUILD_MAP_CFG, CFG_ASSETS_KEY, CFG_GENOMES_KEY

        rgc = RefGenConf.from_yaml_file(map_build_config["config_path"])
        data_dir = rgc.data_dir

        # Find the map config
        map_cfg_path = None
        for genome_dir in os.listdir(data_dir):
            candidate = os.path.join(
                data_dir, genome_dir, "fasta", "default",
                "_refgenie_build", BUILD_MAP_CFG,
            )
            if os.path.exists(candidate):
                map_cfg_path = candidate
                break

        assert map_cfg_path is not None, "Map config not found"

        map_rgc = RefGenConf.from_yaml_file(map_cfg_path)
        assert CFG_GENOMES_KEY in map_rgc
        assert len(map_rgc[CFG_GENOMES_KEY]) > 0

        # Verify fasta asset exists in the map config
        genome_digest = list(map_rgc[CFG_GENOMES_KEY].keys())[0]
        assert "fasta" in map_rgc[CFG_GENOMES_KEY][genome_digest][CFG_ASSETS_KEY]


# =============================================================================
# Test group 2: Reduce step tests
# =============================================================================


class TestReduceBuild:
    """Tests for the --reduce build step."""

    def test_reduce_merges_map_config(self, reduced_config):
        """After reduce, the master config contains the genome and fasta asset."""
        from refgenconf import RefGenConf
        from refgenconf.const import CFG_ASSETS_KEY, CFG_GENOMES_KEY

        rgc = RefGenConf.from_yaml_file(reduced_config["config_path"])
        genomes = rgc.get(CFG_GENOMES_KEY)
        assert genomes is not None and len(genomes) > 0, (
            f"Master config should have genomes after reduce, but got: {genomes}"
        )

        # Verify fasta asset exists
        genome_digest = list(genomes.keys())[0]
        assert "fasta" in genomes[genome_digest][CFG_ASSETS_KEY]

    def test_reduce_preserves_map_config_when_flag_set(self, reduced_config):
        """After reduce with --preserve-map-configs, map config still exists."""
        from refgenconf import RefGenConf
        from refgenconf.const import BUILD_MAP_CFG

        rgc = RefGenConf.from_yaml_file(reduced_config["original_map_config_path"])
        data_dir = rgc.data_dir

        # The map config should still exist since we used --preserve-map-configs
        found = False
        if os.path.exists(data_dir):
            for genome_dir in os.listdir(data_dir):
                map_cfg_path = os.path.join(
                    data_dir, genome_dir, "fasta", "default",
                    "_refgenie_build", BUILD_MAP_CFG,
                )
                if os.path.exists(map_cfg_path):
                    found = True
                    break

        assert found, (
            f"Map config should still exist with --preserve-map-configs"
        )

    def test_reduce_seek_returns_valid_path(self, reduced_config):
        """After reduce, seek returns a valid path to the FASTA file."""
        from refgenconf import RefGenConf
        from refgenconf.const import CFG_GENOMES_KEY

        rgc = RefGenConf.from_yaml_file(reduced_config["config_path"])
        genomes = rgc.get(CFG_GENOMES_KEY)
        assert genomes is not None and len(genomes) > 0

        genome_digest = list(genomes.keys())[0]
        # Use the digest directly since alias may not be available
        fasta_path = rgc.seek(genome_digest, "fasta")
        assert os.path.isfile(fasta_path), (
            f"seek returned nonexistent path: {fasta_path}"
        )


# =============================================================================
# Test group 3: Reduce without --preserve-map-configs
# =============================================================================


class TestReduceWithoutPreserve:
    """Test that reduce without --preserve-map-configs removes map configs."""

    def test_reduce_removes_map_config(self, tmp_path):
        """After reduce without --preserve-map-configs, map config is deleted."""
        config_path = str(tmp_path / "genome_config.yaml")

        # Init
        result = subprocess.run(
            [
                sys.executable, "-m", "refgenie", "init",
                "-c", config_path,
            ],
            capture_output=True, text=True,
        )
        assert result.returncode == 0

        # Map build
        result = subprocess.run(
            [
                sys.executable, "-m", "refgenie", "build",
                "-c", config_path,
                f"{GENOME_ALIAS}/fasta",
                "--map",
                "--files", f"fasta={T7_FASTA_GZ}",
                "--recipe", RECIPE_PARENT,
            ],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, (
            f"Map build failed: {result.stdout}\n{result.stderr}"
        )

        # Verify map config exists before reduce
        from refgenconf import RefGenConf
        from refgenconf.const import BUILD_MAP_CFG

        rgc = RefGenConf.from_yaml_file(config_path)
        data_dir = rgc.data_dir

        map_configs_before = []
        for genome_dir in os.listdir(data_dir):
            candidate = os.path.join(
                data_dir, genome_dir, "fasta", "default",
                "_refgenie_build", BUILD_MAP_CFG,
            )
            if os.path.exists(candidate):
                map_configs_before.append(candidate)

        assert len(map_configs_before) > 0, "Map config not found before reduce"

        # Reduce WITHOUT --preserve-map-configs
        result = subprocess.run(
            [
                sys.executable, "-m", "refgenie", "build",
                "--reduce",
                "-c", config_path,
            ],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, (
            f"Reduce failed: {result.stdout}\n{result.stderr}"
        )

        # Verify map config was removed
        for path in map_configs_before:
            assert not os.path.exists(path), (
                f"Map config should have been removed: {path}"
            )
