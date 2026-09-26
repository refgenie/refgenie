"""
Tests for refgenie.managers.recipe: recipe CRUD, validation, and the
`inherent` declaration round-trip.

Also here: multi-version support in asset classes and recipes. Snakefile
generation is tested in tests/integrations/test_snakemake.py.
"""

from pathlib import Path

import yaml
import pytest
from sqlmodel import Session, select

from refgenie.db.tables import AssetClass, Recipe, RecipeAssetClassesInputs
from refgenie.exceptions import (
    ConfigError,
    MissingAssetClassError,
    MissingRecipeError,
    RecipeExistsError,
)
from refgenie.utils.versioning import parse_name_version, validate_semver


class TestRecipeQuery:
    """list_all / get on registered and missing recipes."""

    def test_minimal_catalog_has_no_recipes(self, refgenie_minimal):
        """A freshly initialized catalog auto-loads no recipes."""
        assert list(refgenie_minimal.recipe.list_all()) == []

    def test_get_registered_recipe(self, refgenie_with_fasta):
        """get returns the fasta recipe registered by the fixture."""
        recipe = refgenie_with_fasta.recipe.get("fasta")
        assert isinstance(recipe, Recipe)
        assert recipe.name == "fasta"

    def test_get_missing_recipe_raises(self, refgenie_with_fasta):
        """get on an unknown name raises MissingRecipeError."""
        with pytest.raises(MissingRecipeError):
            refgenie_with_fasta.recipe.get("nonexistent_recipe_12345")

    @pytest.mark.parametrize("bad_name", [None, ""])
    def test_get_with_empty_name_raises(self, refgenie_with_fasta, bad_name):
        """get with None or an empty name resolves nothing and raises."""
        with pytest.raises(MissingRecipeError):
            refgenie_with_fasta.recipe.get(bad_name)


class TestRecipeAddRemove:
    """Add and remove a recipe, and reject a duplicate add."""

    def test_add_then_remove_round_trip(self, refgenie_with_fasta, fixtures_path):
        """Adding a recipe makes it gettable; removing it makes get raise again."""
        asset_class_file = fixtures_path / "bowtie2_index_asset_class.yaml"
        recipe_file = fixtures_path / "bowtie2_index_asset_recipe.yaml"
        refgenie_with_fasta.asset_class.add(asset_class_file)

        name = yaml.safe_load(recipe_file.read_text())["name"]
        before = len(list(refgenie_with_fasta.recipe.list_all()))

        added = refgenie_with_fasta.recipe.add(recipe_file)
        assert added.name == name
        assert len(list(refgenie_with_fasta.recipe.list_all())) == before + 1
        assert refgenie_with_fasta.recipe.get(name).name == name

        refgenie_with_fasta.recipe.remove(name)
        assert len(list(refgenie_with_fasta.recipe.list_all())) == before
        with pytest.raises(MissingRecipeError):
            refgenie_with_fasta.recipe.get(name)

    def test_add_duplicate_raises(self, refgenie_with_fasta, fixtures_path):
        """Re-adding an already registered recipe raises RecipeExistsError."""
        with pytest.raises(RecipeExistsError):
            refgenie_with_fasta.recipe.add(fixtures_path / "fasta_asset_recipe.yaml")


class TestRecipeValidation:
    """Bad recipe sources and removals surface the right errors."""

    def test_add_nonexistent_file_raises(self, refgenie_minimal, fixtures_path):
        """Adding a path that does not exist raises FileNotFoundError."""
        with pytest.raises(FileNotFoundError):
            refgenie_minimal.recipe.add(fixtures_path / "nonexistent_recipe.yaml")

    def test_add_directory_raises(self, refgenie_minimal, fixtures_path):
        """Adding a directory instead of a file raises IsADirectoryError."""
        with pytest.raises(IsADirectoryError):
            refgenie_minimal.recipe.add(fixtures_path)

    def test_remove_missing_recipe_raises(self, refgenie_minimal):
        """Removing a recipe that does not exist raises MissingRecipeError.

        remove() resolves the recipe (to check for assets built from it) before
        deleting, so a missing name surfaces as MissingRecipeError from get().
        """
        with pytest.raises(MissingRecipeError):
            refgenie_minimal.recipe.remove("nonexistent_recipe_12345")


class TestRecipeInherent:
    """The `inherent` declaration must survive YAML -> database."""

    def test_recipe_inherent_round_trips_from_yaml(
        self, refgenie_with_fasta, fixtures_path, tmp_path
    ):
        """
        A recipe's `inherent` declaration must survive YAML -> database, so
        the digest computed at add-time reflects what the recipe author declared.
        """
        asset_class_file = fixtures_path / "bowtie2_index_asset_class.yaml"
        recipe_file = fixtures_path / "bowtie2_index_asset_recipe.yaml"
        refgenie_with_fasta.asset_class.add(asset_class_file)

        declared = ["*", "!*.log", "!scratch/**"]
        source = yaml.safe_load(recipe_file.read_text())
        source["inherent"] = declared
        declared_file = tmp_path / "bowtie2_with_inherent.yaml"
        declared_file.write_text(yaml.safe_dump(source))

        refgenie_with_fasta.recipe.add(declared_file)
        assert refgenie_with_fasta.recipe.get(source["name"]).inherent == declared

    def test_recipe_without_inherent_is_none(self, refgenie_with_fasta):
        """
        Omitting the declaration must store NULL, not an empty list — None means
        "include everything" while [] would mean "include nothing".
        """
        assert refgenie_with_fasta.recipe.get("fasta").inherent is None


# ---------------------------------------------------------------------------
# Multi-version support in asset classes and recipes
# ---------------------------------------------------------------------------


# --- Version utility tests ---


class TestParseNameVersion:
    @pytest.mark.parametrize(
        "raw, expected",
        [
            ("fasta", ("fasta", None)),
            ("fasta:0.1.0", ("fasta", "0.1.0")),
            ("some:name:1.0.0", ("some:name", "1.0.0")),  # rsplit on last colon
            ("fasta:", ("fasta", "")),
        ],
    )
    def test_parse_name_version(self, raw, expected):
        assert parse_name_version(raw) == expected


class TestValidateSemver:
    @pytest.mark.parametrize(
        "version",
        [
            "0.1.0",
            "1.2.3",
            "10.20.30",
            "1.0.0-alpha",
            "1.0.0-alpha.1",
            "1.0.0-0.3.7",
            "1.0.0+build.1",
            "1.0.0-alpha+001",
        ],
    )
    def test_valid(self, version):
        assert validate_semver(version) is True

    @pytest.mark.parametrize("version", ["foo", "1.2", "1", "", "v1.2.3"])
    def test_invalid(self, version):
        assert validate_semver(version) is False


# --- Asset class multi-version tests ---


def _add_asset_class_directly(session, name, version, description=None):
    """Helper to add an asset class directly to the DB without going through the manager."""
    ac = AssetClass(
        name=name,
        version=version,
        description=description or f"{name} v{version}",
        serving_modes=["file"],
    )
    session.add(ac)
    session.commit()
    session.refresh(ac)
    return ac


def _add_recipe_directly(session, name, version, output_asset_class_id):
    """Helper to add a recipe directly to the DB."""
    recipe = Recipe(
        name=name,
        version=version,
        description=f"{name} v{version}",
        output_asset_class_id=output_asset_class_id,
        command_templates=["echo test"],
        default_asset="default",
    )
    session.add(recipe)
    session.commit()
    session.refresh(recipe)
    return recipe


def _add_asset_class_versions(session, name, versions):
    """Add every version of an asset class named ``name``."""
    for version in versions:
        _add_asset_class_directly(session, name, version)


def _add_recipe_versions(session, name, versions):
    """Add every version of a recipe named ``name`` sharing one output asset class."""
    ac = _add_asset_class_directly(session, f"{name}_ac", "0.1.0")
    for version in versions:
        _add_recipe_directly(session, name, version, ac.id)


class TestAssetClassMultiVersion:
    """Version resolution for the asset_class and recipe managers.

    The resolver body is identical for both managers, so the shared cases are
    parametrized over (manager, version-adder). The asset-class-only cases below
    stay un-parametrized.
    """

    @pytest.mark.parametrize(
        "manager_attr, add_versions",
        [
            ("asset_class", _add_asset_class_versions),
            ("recipe", _add_recipe_versions),
        ],
    )
    def test_get_returns_latest_when_no_version(self, refgenie_minimal, manager_attr, add_versions):
        """get() without version returns the latest by semver."""
        with Session(refgenie_minimal.database_engine) as session:
            add_versions(session, "multitest", ["0.1.0", "0.3.0", "0.2.0"])

        manager = getattr(refgenie_minimal, manager_attr)
        assert manager.get("multitest").version == "0.3.0"

    @pytest.mark.parametrize(
        "manager_attr, add_versions",
        [
            ("asset_class", _add_asset_class_versions),
            ("recipe", _add_recipe_versions),
        ],
    )
    def test_get_returns_specific_version(self, refgenie_minimal, manager_attr, add_versions):
        """get(name, version) returns the specific version."""
        with Session(refgenie_minimal.database_engine) as session:
            add_versions(session, "specific", ["0.1.0", "0.2.0"])

        manager = getattr(refgenie_minimal, manager_attr)
        assert manager.get("specific", "0.1.0").version == "0.1.0"

    @pytest.mark.parametrize(
        "manager_attr, add_versions",
        [
            ("asset_class", _add_asset_class_versions),
            ("recipe", _add_recipe_versions),
        ],
    )
    def test_remove_without_version_errors_when_multiple(
        self, refgenie_minimal, manager_attr, add_versions
    ):
        """remove() without version raises ConfigError when multiple versions exist."""
        with Session(refgenie_minimal.database_engine) as session:
            add_versions(session, "rmtest", ["0.1.0", "0.2.0"])

        manager = getattr(refgenie_minimal, manager_attr)
        with pytest.raises(ConfigError, match="Multiple versions"):
            manager.remove("rmtest")

    @pytest.mark.parametrize(
        "manager_attr, add_versions",
        [
            ("asset_class", _add_asset_class_versions),
            ("recipe", _add_recipe_versions),
        ],
    )
    def test_remove_with_version(self, refgenie_minimal, manager_attr, add_versions):
        """remove(name, version) removes the specific version."""
        with Session(refgenie_minimal.database_engine) as session:
            add_versions(session, "rm_specific", ["0.1.0", "0.2.0"])

        manager = getattr(refgenie_minimal, manager_attr)
        manager.remove("rm_specific", "0.1.0")
        assert not manager.exists("rm_specific", "0.1.0")
        assert manager.exists("rm_specific", "0.2.0")

    def test_get_missing_version_raises(self, refgenie_minimal):
        """get(name, nonexistent_version) raises MissingAssetClassError."""
        with Session(refgenie_minimal.database_engine) as session:
            _add_asset_class_directly(session, "missing_v", "0.1.0")

        with pytest.raises(MissingAssetClassError):
            refgenie_minimal.asset_class.get("missing_v", "9.9.9")

    def test_remove_without_version_succeeds_for_single(self, refgenie_minimal):
        """remove() without version works when only one version exists."""
        with Session(refgenie_minimal.database_engine) as session:
            _add_asset_class_directly(session, "single_rm", "0.1.0")

        refgenie_minimal.asset_class.remove("single_rm")
        assert not refgenie_minimal.asset_class.exists("single_rm", "0.1.0")


class TestRecipeAddWithMultipleAssetClassVersions:
    def test_recipe_add_uses_latest_asset_class(self, refgenie_with_fasta, fixtures_path):
        """recipe add resolves to latest asset class version when multiple exist."""
        # fasta AC v0.1.0 already exists from fixture
        # Add a second version directly
        with Session(refgenie_with_fasta.database_engine) as session:
            _add_asset_class_directly(session, "fasta", "0.2.0")

        # Adding a recipe that references "fasta" should resolve to v0.2.0.
        # The bowtie2 recipe has fasta as an input asset class.
        refgenie_with_fasta.asset_class.add(fixtures_path / "bowtie2_index_asset_class.yaml")
        # Must resolve to one version, not raise MultipleResultsFound
        refgenie_with_fasta.recipe.add(fixtures_path / "bowtie2_index_asset_recipe.yaml")
        recipe = refgenie_with_fasta.recipe.get("bowtie2_index")
        assert recipe is not None


class TestSemverValidation:
    @pytest.mark.parametrize(
        "fixture_name, manager_attr, spec",
        [
            (
                "refgenie_minimal",
                "asset_class",
                {
                    "name": "bad_version_test",
                    "version": "foo",
                    "description": "Test",
                    "serving_modes": ["file"],
                    "seek_keys": {},
                },
            ),
            (
                "refgenie_with_fasta",
                "recipe",
                {
                    "name": "bad_version_recipe",
                    "version": "1.2",
                    "description": "Test",
                    "output_asset_class": "fasta",
                    "command_templates": ["echo test"],
                    "default_asset": "default",
                    "input_params": None,
                    "input_files": None,
                    "input_assets": None,
                    "docker_image": None,
                },
            ),
        ],
    )
    def test_add_rejects_invalid_semver(self, request, fixture_name, manager_attr, spec, tmp_path):
        """Adding a recipe/asset class with an invalid version is rejected."""
        import yaml

        r = request.getfixturevalue(fixture_name)
        bad = tmp_path / "bad_version.yaml"
        bad.write_text(yaml.dump(spec))
        with pytest.raises(ValueError, match="Invalid version"):
            getattr(r, manager_attr).add(bad)


# --- YAML helpers for the on-the-fly asset-class / recipe definitions -------


def _write_yaml(path: Path, data: dict) -> Path:
    """Write a YAML file and return its path."""
    path.write_text(yaml.dump(data, default_flow_style=False))
    return path


def _make_asset_class_yaml(name: str, version: str, **kwargs) -> dict:
    """Create an asset class YAML dict."""
    return {
        "name": name,
        "version": version,
        "description": kwargs.get("description", f"{name} v{version}"),
        "serving_modes": kwargs.get("serving_modes", ["file"]),
        "seek_keys": kwargs.get(
            "seek_keys",
            {
                "main": {
                    "value": "{genome}.dat",
                    "description": "Main file",
                    "type": "file",
                }
            },
        ),
    }


def _make_recipe_yaml(
    name: str,
    version: str,
    output_asset_class: str,
    input_assets: dict | None = None,
) -> dict:
    """Create a recipe YAML dict."""
    return {
        "name": name,
        "version": version,
        "description": f"{name} v{version}",
        "output_asset_class": output_asset_class,
        "command_templates": ["echo test"],
        "default_asset": "default",
        "input_params": None,
        "input_files": None,
        "input_assets": input_assets,
        "docker_image": None,
    }


# =============================================================================
# Asset Class Versioning (integration-only scenarios)
# =============================================================================


class TestAssetClassVersioning:
    """Asset class version scenarios not covered by unit tests."""

    def test_get_latest_with_many_versions(self, refgenie_minimal, tmp_path):
        """Semver ordering picks the correct latest across many versions.

        Exercises numeric vs lexicographic ordering (0.9.0 vs 0.10.0).
        """
        rg = refgenie_minimal

        versions = ["0.1.0", "0.9.0", "0.10.0", "1.0.0", "1.0.1"]
        for v in versions:
            f = _write_yaml(tmp_path / f"ac_{v}.yaml", _make_asset_class_yaml("semver_test", v))
            rg.asset_class.add(f)

        ac = rg.asset_class.get("semver_test")
        assert ac.version == "1.0.1"

    def test_overwrite_existing_version(self, refgenie_minimal, tmp_path):
        """add() with exists_overwrite=True replaces an existing version."""
        rg = refgenie_minimal

        v1 = _write_yaml(
            tmp_path / "ow_v1.yaml",
            _make_asset_class_yaml("overwrite_ac", "0.1.0", description="original"),
        )
        rg.asset_class.add(v1)

        v1_updated = _write_yaml(
            tmp_path / "ow_v1_new.yaml",
            _make_asset_class_yaml("overwrite_ac", "0.1.0", description="updated"),
        )
        rg.asset_class.add(v1_updated, exists_overwrite=True)

        ac = rg.asset_class.get("overwrite_ac", "0.1.0")
        assert ac.description == "updated"


# =============================================================================
# Cross-cutting: Recipe + multi-version asset class
# =============================================================================


class TestRecipeWithMultiVersionAssetClass:
    """Tests that recipe add works when the referenced asset class has multiple versions."""

    def test_recipe_add_with_multiple_output_ac_versions(self, refgenie_minimal, tmp_path):
        """recipe add resolves to latest output asset class version."""
        rg = refgenie_minimal

        # Add two versions of the same asset class
        for v in ["0.1.0", "0.2.0"]:
            f = _write_yaml(
                tmp_path / f"oac_{v}.yaml",
                _make_asset_class_yaml("multi_oac", v),
            )
            rg.asset_class.add(f)

        # Recipe references "multi_oac" without version — should resolve to 0.2.0
        recipe_yaml = _write_yaml(
            tmp_path / "recipe.yaml",
            _make_recipe_yaml("multi_oac_recipe", "0.1.0", "multi_oac"),
        )
        rg.recipe.add(recipe_yaml)

        recipe = rg.recipe.get("multi_oac_recipe")
        assert recipe.output_asset_class.version == "0.2.0"

    def test_recipe_add_with_multiple_input_ac_versions(self, refgenie_minimal, tmp_path):
        """recipe add resolves input asset classes to latest version."""
        rg = refgenie_minimal

        # Output asset class
        out_ac = _write_yaml(
            tmp_path / "out_ac.yaml",
            _make_asset_class_yaml("derived_ac", "0.1.0"),
        )
        rg.asset_class.add(out_ac)

        # Input asset class with two versions
        for v in ["0.1.0", "0.3.0"]:
            f = _write_yaml(
                tmp_path / f"input_ac_{v}.yaml",
                _make_asset_class_yaml("input_ac", v),
            )
            rg.asset_class.add(f)

        # Recipe with input_assets referencing "input_ac" (no version)
        recipe_yaml = _write_yaml(
            tmp_path / "derived_recipe.yaml",
            _make_recipe_yaml(
                "derived_recipe",
                "0.1.0",
                "derived_ac",
                input_assets={
                    "source": {
                        "asset_class": "input_ac",
                        "description": "Source asset",
                        "default": "default",
                    }
                },
            ),
        )
        # Must resolve to one version, not raise MultipleResultsFound
        rg.recipe.add(recipe_yaml)
        recipe = rg.recipe.get("derived_recipe")
        assert recipe is not None

        # The unpinned input asset class must resolve to the LATEST stored
        # version (0.3.0), not the first one added (0.1.0).
        with Session(rg.database_engine) as session:
            link = session.exec(
                select(RecipeAssetClassesInputs).where(
                    RecipeAssetClassesInputs.recipe_id == recipe.id
                )
            ).one()
            resolved = session.get(AssetClass, link.asset_class_id)
        assert resolved.name == "input_ac"
        assert resolved.version == "0.3.0"

    def test_recipe_add_with_pinned_output_ac_version(self, refgenie_minimal, tmp_path):
        """recipe add with 'name:version' syntax pins to a specific AC version."""
        rg = refgenie_minimal

        # Add two versions of the same asset class
        for v in ["0.1.0", "0.2.0"]:
            f = _write_yaml(
                tmp_path / f"pinned_ac_{v}.yaml",
                _make_asset_class_yaml("pinned_ac", v),
            )
            rg.asset_class.add(f)

        # Recipe explicitly pins to v0.1.0
        recipe_yaml = _write_yaml(
            tmp_path / "pinned_recipe.yaml",
            _make_recipe_yaml("pinned_recipe", "0.1.0", "pinned_ac:0.1.0"),
        )
        rg.recipe.add(recipe_yaml)

        recipe = rg.recipe.get("pinned_recipe")
        assert recipe.output_asset_class.version == "0.1.0"
