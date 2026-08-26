"""Tests for asset class and recipe registration."""

from tests.e2e.cli_compat.conftest import assert_exit_error, assert_exit_ok, assert_in_output


class TestAssetClassRecipe:
    def test_add_asset_class(self, initialized_runner, fasta_asset_class_yaml):
        """Add fasta asset class and verify it appears in list."""
        result = initialized_runner.asset_class_add(fasta_asset_class_yaml)
        assert_exit_ok(result)
        result = initialized_runner.asset_class_list()
        assert_exit_ok(result)
        assert_in_output(result, "fasta")

    def test_add_asset_class_duplicate(self, initialized_runner, fasta_asset_class_yaml):
        """Re-adding an existing asset class with force succeeds (overwrite).

        The runner passes -f in Python mode ("already exists. Overwriting");
        the Rust CLI silently re-adds.
        """
        result = initialized_runner.asset_class_add(fasta_asset_class_yaml)
        assert_exit_ok(result)
        result = initialized_runner.asset_class_add(fasta_asset_class_yaml)
        assert_exit_ok(result)
        # The class is still listed after the overwrite
        result = initialized_runner.asset_class_list()
        assert_exit_ok(result)
        assert_in_output(result, "fasta")

    def test_add_recipe(self, initialized_runner, fasta_asset_class_yaml, fasta_recipe_yaml):
        """Add fasta asset class, then recipe. Recipe should appear in list."""
        result = initialized_runner.asset_class_add(fasta_asset_class_yaml)
        assert_exit_ok(result)
        result = initialized_runner.recipe_add(fasta_recipe_yaml)
        assert_exit_ok(result)
        result = initialized_runner.recipe_list()
        assert_exit_ok(result)
        assert_in_output(result, "fasta")

    def test_add_recipe_without_asset_class(self, initialized_runner, fasta_recipe_yaml):
        """Adding a recipe whose output_asset_class does not exist must error.

        Both CLIs exit 1 with a "not found" message (Python: "Asset class
        'fasta' not found."; Rust: "Asset class not found: fasta v*").
        """
        result = initialized_runner.recipe_add(fasta_recipe_yaml)
        assert_exit_error(result)
        assert "not found" in (result.stdout + result.stderr).lower()

    def test_list_asset_classes_empty(self, initialized_runner):
        """On a fresh runner, asset class list should return exit 0."""
        result = initialized_runner.asset_class_list()
        assert_exit_ok(result)

    def test_list_recipes_empty(self, initialized_runner):
        """On a fresh runner, recipe list should return exit 0."""
        result = initialized_runner.recipe_list()
        assert_exit_ok(result)

    def test_add_multiple_asset_classes(
        self,
        initialized_runner,
        fasta_asset_class_yaml,
        bowtie2_asset_class_yaml,
        bwa_asset_class_yaml,
    ):
        """Add fasta, bowtie2_index, and bwa_index asset classes. All should appear in list."""
        for yaml_path in [fasta_asset_class_yaml, bowtie2_asset_class_yaml, bwa_asset_class_yaml]:
            result = initialized_runner.asset_class_add(yaml_path)
            assert_exit_ok(result)
        result = initialized_runner.asset_class_list()
        assert_exit_ok(result)
        assert_in_output(result, "fasta")
        assert_in_output(result, "bowtie2")
        assert_in_output(result, "bwa")

    def test_add_bowtie2_recipe(
        self,
        initialized_runner,
        fasta_asset_class_yaml,
        bowtie2_asset_class_yaml,
        bowtie2_recipe_yaml,
    ):
        """Add bowtie2_index asset class and recipe. Recipe should appear in list."""
        # bowtie2 recipe depends on fasta asset class
        result = initialized_runner.asset_class_add(fasta_asset_class_yaml)
        assert_exit_ok(result)
        result = initialized_runner.asset_class_add(bowtie2_asset_class_yaml)
        assert_exit_ok(result)
        result = initialized_runner.recipe_add(bowtie2_recipe_yaml)
        assert_exit_ok(result)
        result = initialized_runner.recipe_list()
        assert_exit_ok(result)
        assert_in_output(result, "bowtie2")
