"""Tests for build and seek commands."""

import re

import pytest

from tests.e2e.cli_compat.conftest import assert_exit_error, assert_exit_ok


requires_build_tools = pytest.mark.requires_build_tools


class TestBuildSeek:
    def test_seek_nonexistent_asset(self, initialized_runner):
        """Seek a registry path that does not exist. Should error."""
        result = initialized_runner.seek("nonexistent_genome/fasta")
        assert_exit_error(result)

    def test_build_requires_recipe(self, runner_with_genome):
        """Attempt to build without loading the recipe. Should error."""
        result = runner_with_genome.build("demo/fasta")
        assert_exit_error(result)

    @requires_build_tools
    def test_build_fasta_and_seek(self, runner_with_fasta_class):
        """Build fasta for a genome, seek the result, and id the built asset."""
        runner = runner_with_fasta_class
        result = runner.build("demo/fasta", recipe_name="fasta")
        assert_exit_ok(result)
        result = runner.seek("demo/fasta.fasta")
        assert_exit_ok(result)
        # Output should be a file path
        path = result.stdout.strip()
        assert len(path) > 0
        # id on the built asset returns its content digest (64 hex chars).
        result = runner.id("demo/fasta")
        assert_exit_ok(result)
        digest = result.stdout.strip()
        assert re.fullmatch(r"[0-9a-f]{64}", digest), f"unexpected id output: {digest!r}"

    @requires_build_tools
    def test_seek_check_flag_validates_path(self, runner_with_fasta_class):
        """Seek with check should validate the file exists on disk."""
        runner = runner_with_fasta_class
        result = runner.build("demo/fasta", recipe_name="fasta")
        assert_exit_ok(result)
        # Seek with check -- file should exist
        result = runner.seek("demo/fasta.fasta", check=True)
        assert_exit_ok(result)
