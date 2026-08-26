"""Tests for database initialization."""

import pytest

from tests.e2e.cli_compat.conftest import RefgenieRunner, assert_exit_ok


class TestInit:
    def test_init_creates_database(self, runner):
        """Init should succeed and create the database file."""
        result = runner.init()
        assert_exit_ok(result)
        # For python mode, the DB config yaml and sqlite file should exist
        if runner.mode == "python":
            db_config = runner.home_path / "refgenie_db_config.yaml"
            assert db_config.exists(), f"DB config not found at {db_config}"

    def test_init_creates_genome_folder(self, runner):
        """Init should ensure the genome folder directory exists."""
        result = runner.init()
        assert_exit_ok(result)
        assert runner.genome_folder.exists()

    def test_init_idempotent(self, runner):
        """Running init twice should not crash or corrupt the database."""
        result1 = runner.init()
        assert_exit_ok(result1)
        result2 = runner.init()
        assert_exit_ok(result2)
        # After both inits, list genomes should still work
        result3 = runner.list_genomes()
        assert_exit_ok(result3)

    def test_init_nonexistent_parent_dir(self, tmp_path, refgenie_binary, refgenie_mode):
        """Init with a genome_folder whose parent does not exist.

        Both CLIs exit 0. The Python CLI additionally creates the genome
        folder (including parents); refgenie-rs records the folder in its
        config but does not create it yet -- the one observed divergence
        between the two CLIs (2026-08).
        """
        deep_path = tmp_path / "a" / "b" / "c" / "genomes"
        home_path = tmp_path / "home"
        home_path.mkdir()

        r = RefgenieRunner(
            mode=refgenie_mode,
            binary=refgenie_binary,
            db_path=home_path / "refgenie.db",
            genome_folder=deep_path,
            stage_folder=tmp_path / "archives",
            home_path=home_path,
        )
        result = r.init()
        assert_exit_ok(result)  # both CLIs exit 0
        if refgenie_mode == "python":
            assert deep_path.exists(), "init should create genome folder parents"
        else:
            pytest.xfail("refgenie-rs records the genome folder but does not create it yet")
