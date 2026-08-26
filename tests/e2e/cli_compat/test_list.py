"""Tests for list commands.

Only tests with unique coverage live here; genome/asset-class/recipe
listing after adds is covered in test_genome.py and test_recipes.py.
"""

from tests.e2e.cli_compat.conftest import assert_exit_ok, assert_in_output


class TestList:
    def test_list_assets_empty(self, initialized_runner):
        """No assets built. list_assets should return exit 0."""
        result = initialized_runner.list_assets()
        assert_exit_ok(result)

    def test_list_genomes_shows_multiple(self, initialized_runner, test_data_path):
        """Initialize two genomes and both should appear in list."""
        runner = initialized_runner
        result = runner.genome_init("genome1", test_data_path / "demo.fa")
        assert_exit_ok(result)
        result = runner.genome_init("genome2", test_data_path / "rCRSd.fa")
        assert_exit_ok(result)
        result = runner.list_genomes()
        assert_exit_ok(result)
        assert_in_output(result, "genome1")
        assert_in_output(result, "genome2")
