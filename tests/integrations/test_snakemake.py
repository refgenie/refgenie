"""Snakefile generation from refgenie recipes (refgenie/integrations/snakemake/)."""

import re

import yaml

from refgenie.integrations.snakemake.generate import (
    get_snakefile_context,
    populate_snakefile_template,
)


class TestGetSnakefileContext:
    def test_actual_fasta_dependency_preserved(
        self,
        refgenie_with_fasta,
        bowtie2_index_asset_class_file_path,
        bowtie2_index_recipe_file_path,
    ):
        """bowtie2_index declares fasta in input_assets — it still has fasta in input_asset_names."""
        refgenie_with_fasta.asset_class.add(bowtie2_index_asset_class_file_path)
        refgenie_with_fasta.recipe.add(bowtie2_index_recipe_file_path)

        ctx = get_snakefile_context(refgenie_with_fasta)
        bt2_spec = next(s for s in ctx["asset_build_rule_specs"] if s.asset_name == "bowtie2_index")
        assert "fasta" in bt2_spec.input_asset_names

    def test_threads_default(self, refgenie_with_fasta):
        """Every spec has threads=4 in input_params."""
        ctx = get_snakefile_context(refgenie_with_fasta)
        for spec in ctx["asset_build_rule_specs"]:
            assert "threads" in spec.input_params
            assert spec.input_params["threads"] == 4

    def test_recipe_coverage(self, refgenie_with_fasta):
        """asset_names matches the full set of recipe names."""
        ctx = get_snakefile_context(refgenie_with_fasta)
        recipe_names = {r.name for r in refgenie_with_fasta.recipe.list_all()}
        assert set(ctx["asset_names"]) == recipe_names

    def test_second_recipe_version_does_not_duplicate_rule(
        self,
        refgenie_with_fasta,
        bowtie2_index_asset_class_file_path,
        bowtie2_index_recipe_file_path,
        tmp_path,
    ):
        """A recipe stored at two versions still yields exactly ONE rule for it.

        Recipes are immutable: editing one imports a NEW version rather than
        overwriting, so a catalog routinely holds several versions of the same
        recipe. `recipe.list_all()` returns every one of them, so the context
        builder must not take names straight off that list -- that emits the
        name twice and renders two identically-named rules, which makes
        snakemake abort the ENTIRE workflow with "The name build_bowtie2_index
        is already used by another rule". One recipe edit would break every
        build, so this is worth pinning.
        """
        refgenie_with_fasta.asset_class.add(bowtie2_index_asset_class_file_path)
        refgenie_with_fasta.recipe.add(bowtie2_index_recipe_file_path)

        # Import the same recipe again under a higher version, as a real edit would.
        recipe_dict = yaml.safe_load(bowtie2_index_recipe_file_path.read_text())
        recipe_dict["version"] = "9.9.9"
        v2_path = tmp_path / "bowtie2_index_v2.yaml"
        v2_path.write_text(yaml.dump(recipe_dict))
        refgenie_with_fasta.recipe.add(v2_path)

        # Both versions are really in the catalog -- otherwise this proves nothing.
        stored = [
            r.version for r in refgenie_with_fasta.recipe.list_all() if r.name == "bowtie2_index"
        ]
        assert len(stored) == 2, f"expected two stored versions, got {stored}"

        ctx = get_snakefile_context(refgenie_with_fasta)
        assert ctx["asset_names"].count("bowtie2_index") == 1
        specs = [s for s in ctx["asset_build_rule_specs"] if s.asset_name == "bowtie2_index"]
        assert len(specs) == 1


class TestPopulateSnakefileTemplate:
    def test_creates_file(self, refgenie_with_fasta, tmp_path):
        """Output file exists and is non-empty."""
        out = tmp_path / "Snakefile"
        populate_snakefile_template(refgenie_with_fasta, out)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_contains_rules(
        self,
        refgenie_with_fasta,
        tmp_path,
        bowtie2_index_asset_class_file_path,
        bowtie2_index_recipe_file_path,
    ):
        """Output contains the refgenie import, rule all, and a build rule for every
        recipe -- including an added bowtie2_index rule wired to its real fasta and
        genome_init dependencies."""
        refgenie_with_fasta.asset_class.add(bowtie2_index_asset_class_file_path)
        refgenie_with_fasta.recipe.add(bowtie2_index_recipe_file_path)

        out = tmp_path / "Snakefile"
        populate_snakefile_template(refgenie_with_fasta, out)
        content = out.read_text()

        assert "from refgenie import Refgenie" in content
        assert "rule all:" in content
        assert "rule build_fasta:" in content

        recipe_names = [r.name for r in refgenie_with_fasta.recipe.list_all()]
        for name in recipe_names:
            assert f"rule build_{name}:" in content

        # The added bowtie2_index rule references its real (recipe-declared) fasta
        # dependency and the genome_init target.
        assert "rule build_bowtie2_index:" in content
        rule_match = re.search(r"rule build_bowtie2_index:.*?(?=rule |$)", content, re.DOTALL)
        assert rule_match is not None
        rule_block = rule_match.group()
        assert "fasta_target" in rule_block
        assert "genome_init_target" in rule_block

    def test_genome_init_runs_cli_command(self, refgenie_with_fasta, tmp_path):
        """Output has a genome_init rule whose shell command calls refgenie1 genome init."""
        out = tmp_path / "Snakefile"
        populate_snakefile_template(refgenie_with_fasta, out)
        content = out.read_text()
        assert "rule genome_init:" in content
        # Find the genome_init rule block
        rule_match = re.search(r"rule genome_init:.*?(?=rule |$)", content, re.DOTALL)
        assert rule_match is not None
        rule_block = rule_match.group()
        assert "refgenie1 genome init" in rule_block
        assert "--fasta" in rule_block
        assert "--name" in rule_block
        assert "touch" in rule_block

    def test_build_rules_depend_on_genome_init(self, refgenie_with_fasta, tmp_path):
        """Every build rule references genome_init_target in its input."""
        out = tmp_path / "Snakefile"
        populate_snakefile_template(refgenie_with_fasta, out)
        content = out.read_text()

        recipe_names = [r.name for r in refgenie_with_fasta.recipe.list_all()]
        for name in recipe_names:
            rule_match = re.search(rf"rule build_{name}:.*?(?=rule |$)", content, re.DOTALL)
            assert rule_match is not None, f"rule build_{name} not found"
            rule_block = rule_match.group()
            assert "genome_init_target" in rule_block, (
                f"rule build_{name} does not reference genome_init_target"
            )

    def test_no_artificial_fasta_in_non_fasta_rules(self, refgenie_with_fasta, tmp_path):
        """Non-fasta build rules without fasta in recipe input_assets do not reference
        fasta_target. The check is guarded against vacuity: a non-fasta, fasta-free
        recipe is registered so the inner assertion actually runs (refgenie_with_fasta
        alone holds only fasta, which is skipped)."""
        # Register a non-fasta recipe that declares NO fasta input_asset.
        nodep_class = {
            "name": "nodep",
            "version": "0.1.0",
            "description": "no-dependency asset",
            "serving_modes": ["file"],
            "seek_keys": {"out": {"value": "{genome}.out", "description": "out", "type": "file"}},
        }
        nodep_recipe = {
            "name": "nodep",
            "version": "0.1.0",
            "output_asset_class": "nodep",
            "description": "builds nodep with no fasta dependency",
            "input_files": None,
            "input_params": None,
            "input_assets": None,
            "docker_image": None,
            "command_templates": ["echo {{values.output_folder}}"],
            "default_asset": "default",
        }
        class_path = tmp_path / "nodep_class.yaml"
        recipe_path = tmp_path / "nodep_recipe.yaml"
        class_path.write_text(yaml.dump(nodep_class))
        recipe_path.write_text(yaml.dump(nodep_recipe))
        refgenie_with_fasta.asset_class.add(class_path)
        refgenie_with_fasta.recipe.add(recipe_path)

        out = tmp_path / "Snakefile"
        populate_snakefile_template(refgenie_with_fasta, out)
        content = out.read_text()

        recipe_names = [r.name for r in refgenie_with_fasta.recipe.list_all()]
        checked = 0
        for name in recipe_names:
            if name == "fasta":
                continue
            recipe = refgenie_with_fasta.recipe.get(name)
            actual_input_assets = list((recipe.input_assets or {}).keys())
            if "fasta" not in actual_input_assets:
                rule_match = re.search(rf"rule build_{name}:.*?(?=rule |$)", content, re.DOTALL)
                if rule_match:
                    # Extract only the input section
                    input_match = re.search(r"input:.*?output:", rule_match.group(), re.DOTALL)
                    if input_match:
                        input_section = input_match.group()
                        # fasta_target should NOT appear as an input dependency
                        assert "fasta_target" not in input_section, (
                            f"rule build_{name} has artificial fasta_target in input"
                        )
                        checked += 1
        assert checked > 0, "vacuous: no non-fasta, fasta-free rule was actually checked"

    def test_custom_template(self, refgenie_with_fasta, tmp_path):
        """Custom Jinja2 template is rendered correctly."""
        template_file = tmp_path / "custom.smk"
        template_file.write_text("{{ asset_names | join(',') }}")
        out = tmp_path / "Snakefile"
        populate_snakefile_template(refgenie_with_fasta, out, snakefile_template_path=template_file)
        content = out.read_text()

        recipe_names = [r.name for r in refgenie_with_fasta.recipe.list_all()]
        assert content == ",".join(recipe_names)
