"""Tests for remove and tag (rename) commands."""

import pytest

from tests.e2e.cli_compat.conftest import assert_exit_error, assert_exit_ok


requires_build_tools = pytest.mark.requires_build_tools


class TestRemoveTag:
    def test_remove_nonexistent(self, initialized_runner):
        """Remove a registry path that does not exist. Should error."""
        result = initialized_runner.remove("nonexistent/fasta")
        assert_exit_error(result)

    def test_tag_nonexistent(self, initialized_runner):
        """Tag a registry path that does not exist. Should error."""
        result = initialized_runner.tag("nonexistent/fasta:default", "v2")
        assert_exit_error(result)

    @requires_build_tools
    def test_tag_renames_asset(self, runner_with_fasta_class):
        """Build an asset, tag it, verify new name works and old is gone."""
        runner = runner_with_fasta_class
        result = runner.build("demo/fasta", recipe_name="fasta")
        assert_exit_ok(result)

        # Rename default tag to v2
        result = runner.tag("demo/fasta:default", "v2")
        assert_exit_ok(result)

        # Seek new name should succeed
        result = runner.seek("demo/fasta.fasta:v2")
        assert_exit_ok(result)

        # Seek old name should fail
        result = runner.seek("demo/fasta.fasta:default")
        assert_exit_error(result)

    @requires_build_tools
    def test_remove_built_asset(self, runner_with_fasta_class):
        """Build an asset, remove it, verify it is gone."""
        runner = runner_with_fasta_class
        result = runner.build("demo/fasta", recipe_name="fasta")
        assert_exit_ok(result)

        result = runner.remove("demo/fasta")
        assert_exit_ok(result)

        # Seek should fail
        result = runner.seek("demo/fasta.fasta")
        assert_exit_error(result)
