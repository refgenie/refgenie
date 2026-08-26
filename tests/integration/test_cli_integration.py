"""
CLI Integration Tests -- All tests that exercise refgenie via subprocess.

This file covers all CLI subprocess-based integration scenarios.

Run via: ./tests/scripts/test-integration.sh
"""

import gzip
import hashlib
import os
import shutil
import tarfile
from pathlib import Path

import pytest

from tests.helpers import RCRSD_FASTA, TESTS_DATA_DIR, samtools_available
from tests.integration.conftest import (
    TEST_DB_URL,
    _write_db_config,
    run_refgenie,
)

# Skip all integration tests unless explicitly enabled
pytestmark = pytest.mark.skipif(
    os.getenv("RUN_INTEGRATION_TESTS") != "true",
    reason="Integration tests disabled. Run ./tests/scripts/test-integration.sh",
)

# Test genome alias (matches GENOME_ALIAS in conftest.py)
GENOME_ALIAS = "rCRSd"


# =============================================================================
# Helper Functions
# =============================================================================


def _read_fasta_content(path: Path) -> str:
    """Read FASTA content, handling gzip if needed."""
    path = Path(path)
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt") as f:
            return f.read()
    with open(path) as f:
        return f.read()


def _normalize_fasta(text: str) -> dict[str, str]:
    """Parse FASTA into {header: sequence} dict, ignoring line wrapping."""
    sequences = {}
    header = None
    seq_parts = []
    for line in text.strip().splitlines():
        if line.startswith(">"):
            if header is not None:
                sequences[header] = "".join(seq_parts)
            header = line
            seq_parts = []
        else:
            seq_parts.append(line.strip())
    if header is not None:
        sequences[header] = "".join(seq_parts)
    return sequences


# =============================================================================
# Scenario Tests
# =============================================================================


@pytest.mark.shared_state
class TestCLIScenarios:
    """All CLI subprocess integration tests consolidated into scenario tests.

    Each scenario exercises a complete workflow with multiple assertions,
    reusing session-scoped fixtures to avoid redundant builds.
    """

    # -- Scenario 1: CLI Build + Seek + Getseq + ID --

    def test_cli_build_seek_getseq_id(self, cli_build_config):
        """Scenario 1: After build, verify seek/getseq/id all work via CLI."""
        env = cli_build_config["env"]

        # -- Build verification --
        assert cli_build_config["genome_digest"] is not None

        # -- Seek --
        result = run_refgenie("seek", f"{GENOME_ALIAS}/fasta:default", env=env)
        seek_path = result.stdout.strip()
        assert os.path.isfile(seek_path), f"FASTA not found at: {seek_path}"
        assert seek_path.endswith(".fa")

        # Verify content matches original (ignoring line wrapping)
        built_seqs = _normalize_fasta(_read_fasta_content(Path(seek_path)))
        original_seqs = _normalize_fasta(_read_fasta_content(RCRSD_FASTA))
        assert built_seqs == original_seqs

        # -- Getseq --
        result = run_refgenie(
            "getseq",
            "-g",
            GENOME_ALIAS,
            "-l",
            "rCRSd:0-10",
            env=env,
        )
        seq = result.stdout.strip()
        assert len(seq) >= 10
        assert all(c in "ACGTacgt\n" for c in seq)

        # -- ID: genome digest --
        result = run_refgenie("id", GENOME_ALIAS, env=env)
        assert result.returncode == 0
        lines = [x.strip() for x in result.stdout.strip().split("\n") if x.strip()]
        genome_digest = lines[-1]
        assert len(genome_digest) == 32, f"Expected 32-char digest, got: {genome_digest!r}"

        # Genome digest matches fixture
        assert genome_digest == cli_build_config["genome_digest"]

        # -- ID: asset digest with tag --
        result = run_refgenie("id", f"{GENOME_ALIAS}/fasta:default", env=env)
        assert result.returncode == 0
        asset_digest = result.stdout.strip()
        assert len(asset_digest) == 64
        assert all(c in "0123456789abcdef" for c in asset_digest)

        # Genome and asset digests should differ
        assert genome_digest != asset_digest

        # -- List verification via library API --
        from sqlmodel import create_engine
        from refgenie import Refgenie

        engine = create_engine(TEST_DB_URL)
        rg = Refgenie(database_engine=engine, suppress_migrations=True)
        rg.init(genome_folder=Path(cli_build_config["genome_folder"]))
        assets = list(rg.asset.list_assets())
        assert len(assets) >= 1

        # -- Build completion flag --
        # get_asset_build_target_template() advertises the flag at the same path
        # build() writes it to, under the builds/ tree. The advertised path must
        # be the REAL flag, not a symlink to one written elsewhere — otherwise
        # the two roots can drift and snakemake raises MissingOutputException.
        template = rg.get_asset_build_target_template("fasta", "default")
        flag = Path(str(template).replace("{genome_name}", GENOME_ALIAS))
        assert flag.is_file(), f"build flag missing: {flag}"
        assert not flag.is_symlink(), f"build flag should be a real file, not a symlink: {flag}"
        assert f"{os.sep}builds{os.sep}" in str(flag), (
            f"build flag should live under builds/: {flag}"
        )
        # Build bookkeeping must not be inside the asset directory.
        asset_dir = Path(cli_build_config["genome_folder"]) / assets[0].path
        assert not any(p.name.endswith(".flag") for p in asset_dir.rglob("*")), (
            f"build artifacts leaked into the asset directory: {asset_dir}"
        )

        # -- Skip-build path must restore a missing flag --
        # Remove the flag, then rebuild the SAME asset. build() takes the
        # "asset already exists -> skip" early return; it must still restore the
        # declared output, else the nightly fails on previously-built assets
        # with MissingOutputException.
        os.unlink(flag)
        assert not os.path.lexists(flag)
        rg.build_asset(
            recipe_name="fasta",
            genome_name=GENOME_ALIAS,
            asset_group_name="fasta",
            asset_name="default",
        )
        assert flag.is_file(), "skip-build path did not restore the flag"
        assert not flag.is_symlink()
        engine.dispose()

    # -- Scenario 2: Alias Lifecycle --

    def test_cli_alias_lifecycle(self, cli_build_config):
        """Scenario 2: Complete alias lifecycle via CLI."""
        env = cli_build_config["env"]
        genome_digest = cli_build_config["genome_digest"]

        # -- Set alias --
        result = run_refgenie(
            "alias",
            "set",
            "-a",
            "test_alias_flow",
            "-d",
            genome_digest,
            env=env,
        )
        assert result.returncode == 0

        # -- Get alias --
        result = run_refgenie("alias", "get", "-a", "test_alias_flow", env=env)
        assert result.returncode == 0
        assert genome_digest in result.stdout

        # -- Remove alias --
        result = run_refgenie("alias", "remove", "-a", "test_alias_flow", env=env)
        assert result.returncode == 0

        # -- Alias resolves for seek --
        result = run_refgenie("seek", f"{GENOME_ALIAS}/fasta:default", env=env)
        assert result.returncode == 0
        assert result.stdout.strip().endswith(".fa")

        # -- List aliases for genome --
        result = run_refgenie("alias", "get", "-g", genome_digest, env=env)
        assert result.returncode == 0
        assert GENOME_ALIAS in result.stdout

    # -- Scenario 3: Recipe + Asset Class + Multi-Asset Build --

    def test_cli_recipe_asset_class_and_build(self, cli_build_config):
        """Scenario 3: Add asset class, recipe, and build dependent asset via CLI."""
        env = cli_build_config["env"]

        # -- Add asset class via CLI --
        asset_class_path = TESTS_DATA_DIR / "test_asset_class.yaml"
        if not asset_class_path.exists():
            pytest.skip("test_asset_class.yaml not found in test data")
        result = run_refgenie(
            "asset-class", "add", "--source", str(asset_class_path), "-f", env=env
        )
        assert result.returncode == 0

        # -- Verify via library API --
        from sqlmodel import create_engine
        from refgenie import Refgenie

        engine = create_engine(TEST_DB_URL)
        rg = Refgenie(database_engine=engine, suppress_migrations=True)
        rg.init(genome_folder=Path(cli_build_config["genome_folder"]))
        assert rg.recipe.get("fasta") is not None
        assert rg.asset_class.get("fasta") is not None
        engine.dispose()

        # -- Build bwa_index if bwa available --
        if samtools_available() and shutil.which("bwa") is not None:
            recipe_path = TESTS_DATA_DIR / "bwa_index_asset_recipe.yaml"
            if recipe_path.exists():
                run_refgenie(
                    "recipe", "add", "--source", str(recipe_path), "-f", env=env, check=False
                )

            bwa_ac_path = TESTS_DATA_DIR / "bwa_index_asset_class.yaml"
            if bwa_ac_path.exists():
                run_refgenie(
                    "asset-class", "add", "--source", str(bwa_ac_path), "-f", env=env, check=False
                )

            result = run_refgenie(
                "build",
                f"{GENOME_ALIAS}/bwa_index:default",
                "--assets",
                f"fasta={GENOME_ALIAS}/fasta:default",
                env=env,
                check=False,
            )
            if result.returncode == 0:
                engine = create_engine(TEST_DB_URL)
                rg = Refgenie(database_engine=engine, suppress_migrations=True)
                rg.init(genome_folder=Path(cli_build_config["genome_folder"]))
                assert rg.asset.exists(
                    asset_group_name="bwa_index", asset_name="default", genome_name=GENOME_ALIAS
                )
                engine.dispose()

    # -- Scenario 4: Archive + Serve + Download --

    def test_cli_archive_serve_download(self, refgenie_serve_subprocess, tmp_path):
        """Scenario 4: Archive list, real server, genome/asset lists, download, checksum.

        Runs against the real `refgenie serve` subprocess (the same server the
        pull scenario uses), so route or schema drift in refgenie/server/
        fails here. Asset listing uses the real /v4/assets?genome_digest=
        endpoint, and checksum verification uses tarball_digest from the real
        /v4/archives listing -- the same metadata the production puller
        verifies downloads against.
        """
        import httpx

        server = refgenie_serve_subprocess

        # -- Archive verification --
        archive_path = server["archive_path"]
        assert archive_path is not None
        assert os.path.exists(archive_path)

        # Verify tarball contents
        with tarfile.open(archive_path, "r:gz") as tar:
            names = tar.getnames()
            fasta_files = [n for n in names if n.endswith(".fa")]
            fai_files = [n for n in names if n.endswith(".fai")]
            assert len(fasta_files) > 0, "No FASTA files in archive"
            assert len(fai_files) > 0, "No FAI files in archive"

        # Verify database record
        assert server["archive_digest"] is not None
        assert len(server["archive_digest"]) == 64

        # -- Archive list via CLI --
        env = server["env"]
        result = run_refgenie("stage", "list", env=env)
        assert result.returncode == 0
        assert len(result.stdout.strip()) > 0

        # -- Server healthcheck --
        response = httpx.get(f"{server['url']}/v4/healthcheck")
        assert response.status_code == 200
        assert response.json()["status"] == "ok"

        # -- Genome list --
        response = httpx.get(f"{server['url']}/v4/genomes")
        assert response.status_code == 200
        genomes = response.json()["items"]
        assert len(genomes) >= 1
        digests = [g["digest"] for g in genomes]
        assert server["genome_digest"] in digests

        # -- Asset list --
        genome_digest = server["genome_digest"]
        response = httpx.get(
            f"{server['url']}/v4/assets", params={"genome_digest": genome_digest}
        )
        assert response.status_code == 200
        assets = response.json()["items"]
        assert len(assets) >= 1
        asset_groups = [a["asset_group_name"] for a in assets]
        assert "fasta" in asset_groups

        # -- Download + checksum --
        archive_digest = server["archive_digest"]
        resp = httpx.get(f"{server['url']}/v4/archives/{archive_digest}/download")
        assert resp.status_code == 200
        downloaded_content = resp.content
        assert len(downloaded_content) > 0

        listing = httpx.get(
            f"{server['url']}/v4/archives", params={"digest": archive_digest}
        )
        assert listing.status_code == 200
        records = listing.json()["items"]
        assert len(records) == 1
        reported_sha256 = records[0]["tarball_digest"]
        computed_sha256 = hashlib.sha256(downloaded_content).hexdigest()
        assert computed_sha256 == reported_sha256

        # Verify matches local file
        local_content = Path(archive_path).read_bytes()
        assert downloaded_content == local_content

        # -- Extract and verify FASTA content --
        archive_file = tmp_path / "archive.tgz"
        archive_file.write_bytes(resp.content)

        with tarfile.open(archive_file, "r:gz") as tar:
            tar.extractall(path=tmp_path)

        fasta_files = list(tmp_path.rglob("*.fa"))
        assert len(fasta_files) > 0

        extracted_seqs = _normalize_fasta(_read_fasta_content(fasta_files[0]))
        original_seqs = _normalize_fasta(_read_fasta_content(RCRSD_FASTA))
        assert extracted_seqs == original_seqs

    # -- Scenario 5: Pull Operations --

    def test_cli_pull_operations(self, refgenie_serve_subprocess, bulk_pull_client_env):
        """Scenario 5: All pull variants in one test."""
        import httpx

        # -- Server healthcheck --
        response = httpx.get(f"{refgenie_serve_subprocess['url']}/v4/healthcheck")
        assert response.status_code == 200

        # -- Pull single genome --
        result = run_refgenie(
            "pull",
            "-g",
            "rCRSd",
            "--all",
            "--force",
            env=bulk_pull_client_env,
            check=False,
        )
        assert result.returncode == 0, f"pull single genome failed: {result.stderr}"

        # Verify asset was pulled
        seek_result = run_refgenie("seek", "rCRSd/fasta:default", env=bulk_pull_client_env)
        assert os.path.isfile(seek_result.stdout.strip())

        # -- Pull multi genome --
        result = run_refgenie(
            "pull",
            "-g",
            "rCRSd,demo",
            "--all",
            "--force",
            env=bulk_pull_client_env,
            check=False,
        )
        assert result.returncode == 0, (
            f"pull multi genome failed (rc={result.returncode}):\n"
            f"STDOUT: {result.stdout[-1000:]}\n"
            f"STDERR: {result.stderr[-1000:]}"
        )

        # Verify both assets were pulled
        for genome in ["rCRSd", "demo"]:
            seek_result = run_refgenie(
                "seek",
                f"{genome}/fasta:default",
                env=bulk_pull_client_env,
                check=False,
            )
            assert seek_result.returncode == 0, (
                f"seek failed for {genome}: {seek_result.stderr[-500:]}"
            )
            assert os.path.isfile(seek_result.stdout.strip()), f"Asset not found for {genome}"

        # -- Pull asset all genomes --
        result = run_refgenie(
            "pull",
            "--all-genomes",
            "--asset",
            "fasta",
            "--force",
            env=bulk_pull_client_env,
            check=False,
        )
        assert result.returncode == 0, f"pull asset all genomes failed: {result.stderr}"

        # -- Pull init mode --
        result = run_refgenie(
            "pull",
            "-g",
            "rCRSd",
            "--init",
            env=bulk_pull_client_env,
            check=False,
        )
        assert result.returncode == 0, f"pull --init failed: {result.stderr}"

        # -- Mirror --
        result = run_refgenie(
            "mirror",
            "--force",
            env=bulk_pull_client_env,
            check=False,
        )
        assert result.returncode == 0, f"mirror failed: {result.stderr}"

        # -- Pull --all-genomes --all should error --
        result = run_refgenie(
            "pull",
            "--all-genomes",
            "--all",
            env=bulk_pull_client_env,
            check=False,
        )
        assert result.returncode != 0, "pull --all-genomes --all should fail (use mirror)"

    # -- Scenario 6: CLI ID Error Cases --

    def test_cli_id_error_cases(self, cli_env, tmp_path):
        """Scenario 6: CLI id command error handling."""
        test_home = tmp_path / "test_home"
        test_home.mkdir()

        db_config_path = test_home / "refgenie_db_config.yaml"
        _write_db_config(db_config_path, TEST_DB_URL)

        genome_folder = test_home / "genomes"
        genome_folder.mkdir()

        test_env = {
            "REFGENIE_HOME_PATH": str(test_home),
            "REFGENIE_DB_CONFIG_PATH": str(db_config_path),
            "REFGENIE_GENOME_FOLDER": str(genome_folder),
            "NO_COLOR": "1",
            "TERM": "dumb",
            "COLUMNS": "500",
        }

        # Initialize refgenie
        run_refgenie("init", "-f", str(genome_folder), env=test_env)

        # id for nonexistent genome should error
        result = run_refgenie("id", "nonexistent_genome", env=test_env, check=False)
        assert result.returncode != 0

    # -- Scenario 9: Multi-Genome + Compare + List --

    def test_cli_multi_genome_compare_list(self, cli_multi_genome_config):
        """Scenario 9: Multi-genome operations and comparison."""
        from sqlmodel import create_engine
        from refgenie import Refgenie

        env = cli_multi_genome_config["env"]

        # -- Verify both genomes via library API --
        engine = create_engine(TEST_DB_URL)
        rg = Refgenie(database_engine=engine, suppress_migrations=True)
        rg.init(genome_folder=Path(cli_multi_genome_config["genome_folder"]))
        genomes = list(rg.genome.list_all())
        assert len(genomes) >= 2

        # -- Genome compare (same genome) --
        genome_digest = cli_multi_genome_config["genome_digest"]
        result = rg.genome.compare(genome_digest, genome_digest)
        assert result is not None
        assert "digests" in result
        assert result["digests"]["a"] == genome_digest
        assert result["digests"]["b"] == genome_digest

        engine.dispose()

        # -- Genome list via CLI (one subprocess for coverage) --
        result = run_refgenie("genome", "list", env=env)
        assert result.returncode == 0
