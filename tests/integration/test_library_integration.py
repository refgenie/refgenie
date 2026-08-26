"""
Library Integration Tests -- All tests using library API (no CLI subprocesses).

This file covers library API scenarios: read-only operations, mutations,
genome comparison, and bowtie2 child asset operations.

Run via: ./tests/scripts/test-integration.sh
"""

import os
import shutil
from pathlib import Path

import pytest

from tests.helpers import register_fasta

from refgenie import Refgenie
from refgenie.models import BuildParams

# Skip all tests unless explicitly enabled
pytestmark = pytest.mark.skipif(
    os.getenv("RUN_INTEGRATION_TESTS") != "true",
    reason="Integration tests disabled. Run ./tests/scripts/test-integration.sh",
)


def bowtie2_available() -> bool:
    """Check if bowtie2-build is available."""
    return shutil.which("bowtie2-build") is not None


requires_bowtie2 = pytest.mark.skipif(
    not bowtie2_available(), reason="bowtie2-build not available in PATH"
)


def docker_image_available(image: str) -> bool:
    """Check that docker is usable and the given image is present locally."""
    if shutil.which("docker") is None:
        return False
    import subprocess

    if subprocess.run(["docker", "info"], capture_output=True).returncode != 0:
        return False
    return (
        subprocess.run(["docker", "image", "inspect", image], capture_output=True).returncode == 0
    )


requires_bwa_docker = pytest.mark.skipif(
    not docker_image_available("databio/refgenie"),
    reason="docker or databio/refgenie image not available",
)


def build_rcrsd_asset(rg, fixtures_path) -> None:
    """Initialize genome and build the standard rCRSd test asset."""
    register_fasta(rg, fixtures_path)
    rg.genome.initialize_genome(
        fasta_file_path=fixtures_path / "rCRSd.fa",
        alias_names=["rCRSd"],
        description="rCRSd genome",
    )
    rg.build_asset(
        recipe_name="fasta",
        genome_name="rCRSd",
        asset_group_name="fasta",
        asset_name="test",
    )


# =============================================================================
# Scenario 7: Library Read-Only Operations
# =============================================================================


@pytest.mark.shared_state
class TestLibraryReadOnlyOperations:
    """Scenario 7: All read-only library operations in a single class.

    Uses the session-scoped shared_rcrsd_instance to build once. Includes
    self-cleaning mutation tests (alias add/remove, rename/rename-back)
    and server endpoint tests that only read data.
    """

    def test_read_only_operations(self, shared_rcrsd_instance):
        """Scenario 7a: All read-only operations on the built asset."""
        rg = shared_rcrsd_instance

        # -- Seek the asset --
        path = rg.asset.seek("rCRSd", "fasta", "test", force_exists=True)
        assert Path(path).exists()

        # -- Seek again with different API form --
        asset_path = rg.asset.seek(
            genome_name="rCRSd",
            asset_group_name="fasta",
            asset_name="test",
            force_exists=True,
        )
        assert Path(asset_path).exists(), f"FASTA asset file not found: {asset_path}"

        # -- Getseq basic (0-based, half-open) --
        seq = rg.getseq("rCRSd", "rCRSd:0-100")
        assert len(seq) == 100

        # -- Getseq short --
        seq10 = rg.getseq(genome_name="rCRSd", locus="rCRSd:0-10")
        assert isinstance(seq10, str)
        assert len(seq10) == 10
        assert seq10.upper() == seq10

        # -- Genome listing and lookup --
        genomes = list(rg.genome.list_all())
        assert len(genomes) >= 1

        digest = rg.alias.resolve("rCRSd")
        genome = rg.genome.get(digest)
        assert genome.digest == digest
        assert len(digest) == 32

        # -- Asset listing --
        assets = list(rg.asset.list_assets())
        assert len(assets) >= 1

        # -- Asset class and recipe checks --
        assert rg.asset_class.get("fasta") is not None
        assert rg.recipe.get("fasta") is not None

    @pytest.mark.parametrize("locus", ["chr1:0-10", "blah"])
    def test_getseq_raise_errors(self, shared_rcrsd_instance, locus):
        """Scenario 7b: getseq raises errors for invalid loci."""
        rg = shared_rcrsd_instance
        with pytest.raises(ValueError):
            rg.getseq(genome_name="rCRSd", locus=locus)

    def test_self_cleaning_asset_rename(self, shared_rcrsd_instance):
        """Scenario 7d: Asset rename and rename-back (self-cleaning)."""
        rg = shared_rcrsd_instance
        genome_digest = rg.alias.resolve("rCRSd")

        # Get the original asset
        original_asset = rg.asset.get(
            genome_digest=genome_digest, asset_group_name="fasta", asset_name="test"
        )
        assert original_asset is not None
        assert original_asset.name == "test"

        # Check the group default before rename (it lives on the AssetName flag)
        original_default_asset = rg.asset.get_default("fasta", genome_digest=genome_digest)

        # Rename the asset
        new_name = "test_renamed"
        renamed_asset = rg.asset.rename(
            genome_digest=genome_digest,
            asset_group_name="fasta",
            asset_name="test",
            new_asset_name=new_name,
        )

        # Verify the asset was renamed
        assert renamed_asset is not None
        assert renamed_asset.name == new_name
        assert renamed_asset.digest == original_asset.digest

        # Verify the original asset name no longer exists
        assert not rg.asset.exists("fasta", "test", genome_name="rCRSd")

        # Verify the new asset name exists
        assert rg.asset.exists("fasta", new_name, genome_name="rCRSd")

        # Renaming the default name carries the default flag with the row.
        if original_default_asset == "test":
            assert rg.asset.get_default("fasta", genome_digest=genome_digest) == new_name

        # Rename back to original name for cleanup
        rg.asset.rename(
            genome_digest=genome_digest,
            asset_group_name="fasta",
            asset_name=new_name,
            new_asset_name="test",
        )

        # Verify it was renamed back
        assert rg.asset.exists("fasta", "test", genome_name="rCRSd")
        assert not rg.asset.exists("fasta", new_name, genome_name="rCRSd")

    def test_server_endpoints_with_built_asset(self, shared_client, shared_rcrsd_instance):
        """Scenario 7e: Server endpoints that read from a built asset.

        Uses shared_client (backed by same SQLite DB as shared_rcrsd_instance).
        """
        # -- Genome --
        digest = shared_rcrsd_instance.alias.resolve("rCRSd")
        response = shared_client.get(f"/v4/genomes/{digest}")
        assert response.status_code == 200

        # -- Summary --
        response = shared_client.get("/v4/summary")
        assert response.status_code == 200
        data = response.json()
        assert data["genomes"] >= 1
        assert data["assets"] >= 1


# =============================================================================
# Scenario 8: Library Mutation Operations
# =============================================================================


class TestLibraryMutationOperations:
    """Scenario 8: Library mutation operations that need fresh DB state.

    Each test gets a fresh database via the clean_database autouse fixture.
    """

    def test_genome_reinit_and_asset_removal(self, test_engine, fixtures_path, tmp_path):
        """Scenario 8: Genome reinit rejection + asset removal."""
        rg = Refgenie(database_engine=test_engine, suppress_migrations=True)
        rg.init(genome_folder=tmp_path / "genomes")
        register_fasta(rg, fixtures_path)

        fasta_path = fixtures_path / "rCRSd.fa"

        # -- Reinit rejection --
        rg.genome.initialize_genome(
            fasta_file_path=fasta_path,
            alias_names=["rCRSd"],
            description="rCRSd genome",
        )

        with pytest.raises(ValueError):
            rg.genome.initialize_genome(
                fasta_file_path=fasta_path,
                alias_names=["genome2"],
                description="rCRSd genome",
            )

        # -- Asset build and removal --
        rg.build_asset(
            recipe_name="fasta",
            genome_name="rCRSd",
            asset_group_name="fasta",
            asset_name="default",
        )

        digest = rg.alias.resolve("rCRSd")

        assert rg.asset.exists(
            genome_name="rCRSd",
            asset_group_name="fasta",
            asset_name="default",
        )

        rg.asset.remove(
            genome_name="rCRSd",
            asset_group_name="fasta",
            asset_name="default",
        )

        # Asset should no longer exist. Removing a genome's last asset cascades:
        # the empty group is removed, and with no groups left the genome (and its
        # aliases) go too -- so the alias no longer resolves and we query by the
        # digest captured before removal.
        assert not rg.asset.exists(
            genome_digest=digest,
            asset_group_name="fasta",
            asset_name="default",
        )

        # Removing the last asset removes the now-empty genome as well.
        assert not rg.genome.exists(digest)

        # Verify list_assets returns empty
        assert not rg.asset.list_assets()


# =============================================================================
# Scenario: Two Genomes Compare
# =============================================================================


@pytest.mark.shared_state
class TestTwoGenomesCompare:
    """Test genome comparison using class-scoped fixture with two genomes.

    Uses its own SQLite database to avoid conflicting with the shared
    session-scoped engine used by read-only tests.
    """

    @pytest.fixture(autouse=True, scope="class")
    def two_genome_instance(self, fixtures_path, tmp_path_factory):
        """Build two genomes for comparison tests with isolated SQLite DB."""
        from sqlmodel import create_engine, SQLModel

        db_dir = tmp_path_factory.mktemp("two_genome_db")
        db_path = db_dir / "two_genome.db"
        engine = create_engine(f"sqlite:///{db_path}")
        SQLModel.metadata.create_all(engine)

        tmp_dir = tmp_path_factory.mktemp("two_genome_compare")
        rg = Refgenie(database_engine=engine, suppress_migrations=True)
        rg.init(genome_folder=tmp_dir / "genomes")
        register_fasta(rg, fixtures_path)

        rg.genome.initialize_genome(
            fasta_file_path=fixtures_path / "rCRSd.fa",
            alias_names=["rCRSd"],
            description="rCRSd mitochondrial reference",
        )
        rg.build_asset(
            recipe_name="fasta",
            genome_name="rCRSd",
            asset_group_name="fasta",
            asset_name="v1",
        )

        rg.genome.initialize_genome(
            fasta_file_path=fixtures_path / "t7.fa",
            alias_names=["t7"],
            description="T7 phage genome",
        )
        rg.build_asset(
            recipe_name="fasta",
            genome_name="t7",
            asset_group_name="fasta",
            asset_name="v1",
        )

        yield rg
        engine.dispose()

    def test_two_genomes_and_compare(self, two_genome_instance):
        """Compare two pre-built genomes."""
        rg = two_genome_instance

        genomes = list(rg.genome.list_all())
        assert len(genomes) == 2

        digest_a = rg.alias.resolve("rCRSd")
        digest_b = rg.alias.resolve("t7")
        result = rg.genome.compare(digest_a, digest_b)

        assert result is not None
        assert isinstance(result, dict)


# =============================================================================
# Scenario: Bowtie2 Operations
# =============================================================================


@requires_bowtie2
@pytest.mark.shared_state
class TestBowtie2Operations:
    """Bowtie2 child asset operations with shared class fixture.

    Builds fasta + bowtie2 once for the class.
    """

    @pytest.fixture(autouse=True, scope="class")
    def bowtie2_setup(self, fixtures_path, tmp_path_factory):
        """Build fasta + bowtie2 once for the class with isolated SQLite DB."""
        from sqlmodel import create_engine, SQLModel

        db_dir = tmp_path_factory.mktemp("bowtie2_db")
        db_path = db_dir / "bowtie2.db"
        engine = create_engine(f"sqlite:///{db_path}")
        SQLModel.metadata.create_all(engine)

        tmp_dir = tmp_path_factory.mktemp("bowtie2_genomes")
        rg = Refgenie(database_engine=engine, suppress_migrations=True)
        rg.init(genome_folder=tmp_dir / "genomes")
        build_rcrsd_asset(rg, fixtures_path)

        rg.asset_class.add(fixtures_path / "bowtie2_index_asset_class.yaml")
        rg.recipe.add(fixtures_path / "bowtie2_index_asset_recipe.yaml")
        rg.build_asset(
            recipe_name="bowtie2_index",
            genome_name="rCRSd",
            asset_group_name="bowtie2_index",
            asset_name="test",
            params=BuildParams(params={"threads": 2}),
        )

        yield rg
        engine.dispose()

    def test_bowtie2_parent_child_operations(self, bowtie2_setup):
        """Bowtie2 child asset: build, verify parent/child, can't remove parent."""
        rg = bowtie2_setup

        # -- Verify bowtie2 asset exists --
        assert rg.asset.exists(
            asset_group_name="bowtie2_index", asset_name="test", genome_name="rCRSd"
        )

        d = rg.alias.resolve("rCRSd")

        # -- Parent/child relationships --
        assert rg.asset.get_parents(
            genome_digest=d, asset_group_name="bowtie2_index", asset_name="test"
        )
        assert rg.asset.get_children(asset_group_name="fasta", asset_name="test", genome_digest=d)

        # -- Can't remove parent with children --
        with pytest.raises(ValueError):
            rg.asset.remove(genome_name="rCRSd", asset_group_name="fasta", asset_name="test")


# =============================================================================
# Scenario: Docker build of a colocation-dependent child asset
# =============================================================================


@requires_bwa_docker
class TestDockerColocationBuild:
    """Build bwa_index (which colocates the parent fasta) inside a docker container.

    Regression test for the docker colocation-mount bug: the builder created a
    *relative* colocation symlink (e.g. bwa_index/default/<digest>.fa ->
    ../../fasta/default/<digest>.fa) but mounted ONLY the child output folder
    into the container, so the symlink target dangled inside the container and
    `bwa index <digest>.fa` failed with "No such file or directory".

    The fix is to also mount the genome folder (which contains both the child
    output and the parent fasta asset) so relative symlinks resolve in-container.
    """

    BWA_INDEX_FILES = (".amb", ".ann", ".bwt", ".pac", ".sa")

    @pytest.fixture(scope="class")
    def bwa_docker_setup(self, fixtures_path, tmp_path_factory):
        """Init genome, build fasta natively, then build bwa_index with docker=True."""
        from sqlmodel import create_engine, SQLModel

        db_dir = tmp_path_factory.mktemp("bwa_docker_db")
        db_path = db_dir / "bwa_docker.db"
        engine = create_engine(f"sqlite:///{db_path}")
        SQLModel.metadata.create_all(engine)

        tmp_dir = tmp_path_factory.mktemp("bwa_docker_genomes")
        rg = Refgenie(database_engine=engine, suppress_migrations=True)
        rg.init(genome_folder=tmp_dir / "genomes")

        # Build the parent fasta asset (native; no docker needed).
        build_rcrsd_asset(rg, fixtures_path)

        # Register the bwa_index asset class + recipe (recipe colocates the fasta).
        rg.asset_class.add(fixtures_path / "bwa_index_asset_class.yaml")
        rg.recipe.add(fixtures_path / "bwa_index_asset_recipe.yaml")

        # Build the colocation-dependent child asset INSIDE docker (the
        # regression path: a dangling colocation symlink).
        asset = rg.build_asset(
            recipe_name="bwa_index",
            genome_name="rCRSd",
            asset_group_name="bwa_index",
            asset_name="test",
            docker=True,
        )

        yield rg, asset
        engine.dispose()

    def test_docker_colocation_build_produces_index(self, bwa_docker_setup):
        """The docker bwa_index build must succeed and record the asset + outputs."""
        rg, asset = bwa_docker_setup

        # Build must not have failed/skipped (None means failure).
        assert asset is not None, (
            "Docker build of colocation-dependent bwa_index returned None "
            "(build failed -- likely dangling colocation symlink in container)."
        )

        # Asset must be recorded in the database.
        assert rg.asset.exists(asset_group_name="bwa_index", asset_name="test", genome_name="rCRSd")

        # The bwa index output files must exist on disk. The data folder keys an
        # asset's output directory by its content digest, not the asset name.
        digest = rg.alias.resolve("rCRSd")
        output_folder = rg.data_folder / digest / "bwa_index" / asset.digest
        fa_base = output_folder / f"{digest}.fa"
        missing = [
            ext
            for ext in self.BWA_INDEX_FILES
            if not (output_folder / f"{digest}.fa{ext}").exists()
        ]
        assert not missing, (
            f"Missing bwa index output files {missing} in {output_folder}. "
            f"Existing: {sorted(p.name for p in output_folder.iterdir())}"
        )

        # The colocation symlink itself should resolve to the parent fasta.
        assert fa_base.exists(), f"Colocation symlink does not resolve: {fa_base}"
