"""
CLI compatibility test fixtures.

Provides RefgenieRunner abstraction that works against either
the Python CLI or the Rust CLI binary. Tests use only subprocess
calls -- no Python library imports from refgenie.

Mode selection is explicit: the suite runs against the Python CLI unless
REFGENIE_MODE=rust is set. To run against the Rust CLI:

    REFGENIE_MODE=rust REFGENIE_BIN=/path/to/refgenie pytest -m e2e tests/e2e/cli_compat/
"""

import os
import shutil

import pytest

from tests.e2e.cli_compat.helpers import RefgenieRunner, assert_exit_ok


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "requires_build_tools: mark test as requiring external build tools (samtools, etc.)",
    )


def pytest_collection_modifyitems(config, items):
    """Skip requires_build_tools tests if samtools is not available (Python mode only).

    In Rust mode, build tools like samtools are not needed because the Rust
    binary (refgenie-build-fasta) does FASTA indexing natively.
    """
    mode = os.environ.get("REFGENIE_MODE", "").lower()
    if mode == "rust":
        return
    if not shutil.which("samtools"):
        skip_marker = pytest.mark.skip(reason="samtools not available in PATH (Python mode)")
        for item in items:
            if "requires_build_tools" in item.keywords:
                item.add_marker(skip_marker)


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def refgenie_binary():
    """Resolve the REFGENIE_BIN path."""
    binary = os.environ.get("REFGENIE_BIN", "refgenie")
    # For python mode we use sys.executable directly in the runner
    return binary


@pytest.fixture(scope="session")
def refgenie_mode():
    """Return "rust" if REFGENIE_MODE=rust, otherwise "python".

    There is no auto-detection: the Rust binary's --version output
    ("refgenie 0.1.0") is indistinguishable from the Python CLI's, so
    rust mode must be requested explicitly (see module docstring).
    """
    return "rust" if os.environ.get("REFGENIE_MODE", "").lower() == "rust" else "python"


@pytest.fixture(scope="session")
def test_data_path(fixtures_path):
    """Alias for the root ``fixtures_path`` fixture.

    The ``test_``-prefixed name is a collection hazard (pytest very nearly
    treats it as a test); it survives only because the cli_compat test modules
    still ask for it by that name.
    """
    return fixtures_path


@pytest.fixture(scope="session")
def fasta_asset_class_yaml(fixtures_path):
    return fixtures_path / "fasta_asset_class.yaml"


@pytest.fixture(scope="session")
def fasta_recipe_yaml(fixtures_path):
    return fixtures_path / "fasta_asset_recipe.yaml"


@pytest.fixture(scope="session")
def bowtie2_asset_class_yaml(fixtures_path):
    return fixtures_path / "bowtie2_index_asset_class.yaml"


@pytest.fixture(scope="session")
def bowtie2_recipe_yaml(fixtures_path):
    return fixtures_path / "bowtie2_index_asset_recipe.yaml"


@pytest.fixture(scope="session")
def bwa_asset_class_yaml(fixtures_path):
    return fixtures_path / "bwa_index_asset_class.yaml"


@pytest.fixture(scope="session")
def demo_fasta(fixtures_path):
    return fixtures_path / "demo.fa"


@pytest.fixture
def runner(tmp_path, refgenie_binary, refgenie_mode):
    """Create a fresh RefgenieRunner with isolated tmp_path database and genome folder."""
    home_path = tmp_path / "refgenie_home"
    home_path.mkdir()
    genome_folder = tmp_path / "genomes"
    genome_folder.mkdir()
    stage_folder = tmp_path / "archives"
    stage_folder.mkdir()
    db_path = home_path / "refgenie.db"

    return RefgenieRunner(
        mode=refgenie_mode,
        binary=refgenie_binary,
        db_path=db_path,
        genome_folder=genome_folder,
        stage_folder=stage_folder,
        home_path=home_path,
    )


@pytest.fixture
def initialized_runner(runner):
    """Runner with init already called."""
    result = runner.init()
    assert_exit_ok(result)
    return runner


@pytest.fixture
def runner_with_genome(initialized_runner, demo_fasta):
    """Runner with init called and a demo genome initialized."""
    result = initialized_runner.genome_init("demo", demo_fasta)
    assert_exit_ok(result)
    return initialized_runner


@pytest.fixture
def runner_with_fasta_class(runner_with_genome, fasta_asset_class_yaml, fasta_recipe_yaml):
    """Runner with fasta asset class and recipe loaded, plus a demo genome."""
    result = runner_with_genome.asset_class_add(fasta_asset_class_yaml)
    assert_exit_ok(result)
    result = runner_with_genome.recipe_add(fasta_recipe_yaml)
    assert_exit_ok(result)
    return runner_with_genome
