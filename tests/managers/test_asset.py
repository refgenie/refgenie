"""
Tests for the asset manager CRUD/registry surface (refgenie.managers.asset).

Write-path correctness (content.add atomicity, write ordering) lives in
test_asset_content.py; removal lives in test_asset_removal.py. This file
covers registry-path parsing/population, the read-only query surface
(exists/seek/get/list/seek-keys/table), the AssetClass registration lifecycle,
and the unit half of seek-key handling (the built-in seek-keys builder, CLI
parsing, and asset-class seek-key persistence).

The component half of asset names lives in test_asset_content.py; the
component half of seek keys (non-path seek keys) lives below in this file.
Remote seek (seekr) lives in test_servers.py.
"""

import json
import shutil
import sys
from io import StringIO
from pathlib import Path
from subprocess import CalledProcessError, TimeoutExpired
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest
from pydantic import ValidationError

from refgenie import Refgenie
from refgenie.cli.commands.curate import handle_add
from refgenie.db.tables import SeekKeyType
from refgenie.exceptions import (
    AssetClassExistsError,
    CustomSeekKeyError,
    MissingAliasError,
    MissingAssetClassError,
    MissingAssetError,
)
from refgenie.managers.asset import AssetManager
from refgenie.managers.build import CUSTOM_SEEK_KEY_TIMEOUT, BuildManager
from refgenie.managers.sources import ServerManager
from tests.helpers import (
    add_asset_from_files,
    fake_digest,
    make_command_values,
    make_engine,
    register_fasta,
)


class TestAssetRegistryPath:
    """Parsing and population of refgenie:// registry paths (unit tier)."""

    @pytest.mark.parametrize(
        "asset_registry_path, result",
        [
            (
                "hg38/fasta:default",
                {
                    "protocol": None,
                    "genome": "hg38",
                    "asset_group": "fasta",
                    "seek_key": None,
                    "asset": "default",
                },
            ),
            (
                "hg38/fasta:custom",
                {
                    "protocol": None,
                    "genome": "hg38",
                    "asset_group": "fasta",
                    "seek_key": None,
                    "asset": "custom",
                },
            ),
            (
                "hg38/fasta.fasta:custom",
                {
                    "protocol": None,
                    "genome": "hg38",
                    "asset_group": "fasta",
                    "seek_key": "fasta",
                    "asset": "custom",
                },
            ),
            (
                "fasta",
                {
                    "protocol": None,
                    "genome": None,
                    "asset_group": "fasta",
                    "seek_key": None,
                    "asset": None,
                },
            ),
            (
                "hg38/fasta",
                {
                    "protocol": None,
                    "genome": "hg38",
                    "asset_group": "fasta",
                    "seek_key": None,
                    "asset": None,
                },
            ),
            (
                "test_proto://hg38/fasta",
                {
                    "protocol": "test_proto",
                    "genome": "hg38",
                    "asset_group": "fasta",
                    "seek_key": None,
                    "asset": None,
                },
            ),
            (
                "refgenie://hg38/fasta.fasta:custom",
                {
                    "protocol": "refgenie",
                    "genome": "hg38",
                    "asset_group": "fasta",
                    "seek_key": "fasta",
                    "asset": "custom",
                },
            ),
        ],
    )
    def test_valid_parsing(self, asset_registry_path, result):
        assert Refgenie.parse_asset_registry_path(asset_registry_path).model_dump() == result

    @pytest.mark.parametrize(
        "asset_registry_path",
        [
            "hg38/fasta:defa%-",
            "hg38/fasta:custo$#",
            "hg38/fasta:custom:tag",
            "hg38/fasta/test:default",
        ],
    )
    def test_invalid_tag_raises(self, asset_registry_path):
        with pytest.raises(ValidationError):
            Refgenie.parse_asset_registry_path(asset_registry_path)

    @pytest.mark.parametrize(
        "asset_registry_path",
        [
            "hg38/fasta.fasta:default",
            "unknown://hg38/fasta:default",
        ],
    )
    def test_populating_requires_refgenie_protocol(self, refgenie_minimal, asset_registry_path):
        """A non-refgenie:// path is returned untouched."""
        with patch.object(AssetManager, "seek_components", return_value=Path("whatever")):
            assert refgenie_minimal.populate(asset_registry_path) == asset_registry_path

    @pytest.mark.parametrize(
        "asset_registry_path, result, value",
        [
            ("refgenie://hg38/fasta.fasta:default", "REPLACED", "REPLACED"),
            ("refgenie://hg38/fasta:default", "REPLACED", "REPLACED"),
            ({"key": "refgenie://hg38/fasta:default"}, {"key": "REPLACED"}, "REPLACED"),
            (
                {"key": "refgenie://hg38/fasta:default", "key2": "refgenie://hg38/fasta:default"},
                {"key": "REPLACED", "key2": "REPLACED"},
                "REPLACED",
            ),
            (["refgenie://hg38/fasta:default"], ["REPLACED"], "REPLACED"),
            (
                ["ABC refgenie://hg38/fasta:default XYZ", "refgenie://hg38/fasta:default"],
                ["ABC REPLACED XYZ", "REPLACED"],
                "REPLACED",
            ),
            (
                {
                    "key1": "refgenie://hg38/fasta:default",
                    "key2": {"key": "refgenie://hg38/fasta:default"},
                },
                {"key1": "REPLACED", "key2": {"key": "REPLACED"}},
                "REPLACED",
            ),
        ],
    )
    def test_populating(self, refgenie_minimal, asset_registry_path, result, value):
        """refgenie:// paths are replaced across str/dict/list/nested/embedded inputs."""
        with patch.object(AssetManager, "seek_components", return_value=str(Path(value))):
            assert refgenie_minimal.populate(asset_registry_path) == result

        with patch.object(ServerManager, "seek_components", return_value=value):
            assert refgenie_minimal.populater(asset_registry_path) == result


class TestAssetQuery:
    """Read-only asset query surface on a built catalog (unit tier)."""

    def test_exists(self, refgenie_session):
        """exists() is True for the built asset and False for each way to miss;
        get() returns the populated asset with its name and a non-empty digest."""
        r = refgenie_session
        digest = r.alias.resolve("rCRSd")
        assert r.asset.exists("fasta", "test", genome_digest=digest)
        assert not r.asset.exists("nonexistent", "test", genome_digest=digest)
        assert not r.asset.exists("fasta", "nonexistent", genome_digest=digest)
        assert not r.asset.exists("fasta", "test", genome_digest=fake_digest("nonexistent"))

        asset = r.asset.get(genome_digest=digest, asset_group_name="fasta", asset_name="test")
        assert asset.name == "test"
        assert isinstance(asset.digest, str) and asset.digest

    def test_seek(self, refgenie_session):
        """seek returns a JSON-serializable str path to a real file; missing raises."""
        r = refgenie_session
        digest = r.alias.resolve("rCRSd")
        path = r.asset.seek(digest, "fasta", "test")
        assert isinstance(path, str)
        json.dumps(path)  # seek values must be JSON-serializable
        assert Path(path).exists()
        # force_exists returns the same real path.
        assert r.asset.seek(digest, "fasta", "test", force_exists=True) == path
        with pytest.raises(MissingAssetError):
            r.asset.seek(digest, "nonexistent", "test", force_exists=True)

    def test_seek_by_digest_matches_registry_path(self, refgenie_session):
        """seek(<digest>) resolves the same as seek_components(<alias>/fasta).

        The genome has a local alias, so a digest falls back to that alias tree
        and yields an identical path.
        """
        r = refgenie_session
        digest = r.alias.resolve("rCRSd")
        assert digest != "rCRSd"
        by_alias = r.asset.seek_components(r.parse_asset_registry_path("rCRSd/fasta:test"))
        by_digest = r.asset.seek(digest, "fasta", "test")
        assert by_digest == by_alias

    def test_registry_path_genome_is_an_alias(self, refgenie_session):
        """A digest in a registry path is not looked up as a digest."""
        r = refgenie_session
        digest = r.alias.resolve("rCRSd")
        with pytest.raises(MissingAliasError):
            r.asset.seek_components(r.parse_asset_registry_path(f"{digest}/fasta:test"))

    def test_list_assets_by_digest(self, refgenie_session):
        """list_assets (behind `list -g` and `list --genome-digest`) filters by digest."""
        r = refgenie_session
        digest = r.alias.resolve("rCRSd")
        assets = list(r.asset.list_assets(genome_digests=[digest]))
        assert assets
        assert "fasta" in {a.asset_group.name for a in assets}

    def test_list_assets_and_groups(self, refgenie_session):
        """list_assets returns the built asset(s), and the fasta group is among them."""
        assets = list(refgenie_session.asset.list_assets())
        assert assets
        assert "fasta" in {a.asset_group.name for a in assets}

    def test_list_seek_keys_values_shape(self, refgenie_session):
        """Bulk list_seek_keys_values returns a nested, str-only, JSON-serializable mapping."""
        result = refgenie_session.asset.list_seek_keys_values()
        assert isinstance(result, dict)
        seek_map = result["rCRSd"]["fasta"]["test"]
        assert isinstance(seek_map, dict) and seek_map
        # All values must be str and JSON-serializable.
        for sk_name, sk_val in seek_map.items():
            assert isinstance(sk_name, str)
            assert isinstance(sk_val, str)
        json.dumps(result)

    def test_list_seek_keys_values_filter_by_asset_group(self, refgenie_session):
        """list_seek_keys_values respects the asset_group_name filter."""
        result = refgenie_session.asset.list_seek_keys_values(asset_group_name="fasta")
        assert "rCRSd" in result
        for _genome_name, groups in result.items():
            assert set(groups.keys()) == {"fasta"}

    def test_assets_table_regression(self, refgenie_session):
        """asset.table() returns a list of rich Tables.

        Regression: the AliasManager.list() -> list_all() rename left a stale
        call in AssetManager, crashing `refgenie list` with AttributeError.
        """
        from rich.table import Table

        tables = refgenie_session.asset.table()
        assert isinstance(tables, list) and tables
        assert all(isinstance(t, Table) for t in tables)


# --- Asset-class registration lifecycle (unit tier) --------------------------


class TestAssetClassSchema:
    """AssetClass registration lifecycle and duplicate contract (unit tier)."""

    def test_add_get_remove_lifecycle(self, refgenie_minimal, fixtures_path):
        """Add increments the count and is gettable; remove reverses both and get then raises."""
        r = refgenie_minimal
        initial = len(list(r.asset_class.list_all()))

        r.asset_class.add(fixtures_path / "bowtie2_index_asset_class.yaml")
        assert len(list(r.asset_class.list_all())) == initial + 1
        assert r.asset_class.get("bowtie2_index").name == "bowtie2_index"

        r.asset_class.remove("bowtie2_index")
        assert len(list(r.asset_class.list_all())) == initial
        with pytest.raises(MissingAssetClassError):
            r.asset_class.get("bowtie2_index")

    def test_duplicate_add_raises(self, refgenie_with_fasta, fixtures_path):
        """Re-adding an already-registered asset class raises AssetClassExistsError."""
        with pytest.raises(AssetClassExistsError):
            refgenie_with_fasta.asset_class.add(fixtures_path / "fasta_asset_class.yaml")


class TestEmptyDefaultAssetNameIsFatal:
    """An unresolvable tool version must raise, never silently become "default" (unit tier).

    Assets are named after the version of the tool that built them, produced by a
    shell pipeline (grep/awk/head) that yields an EMPTY string with a zero exit
    status when the tool fails to report a version. `or "default"` turned that into
    a published, mis-named asset baked into the public S3 key, and made callers'
    "refuse to fall back" guards dead code by handing them a truthy string.
    """

    @staticmethod
    def _namespaces(version=""):
        return make_command_values(
            custom_seek_keys={"version": version},
            asset_group_name="hisat2_index",
            genome_digest="",
            genome_folder=Path(""),
        )

    @pytest.mark.parametrize("version", ["", "   \n"])
    def test_empty_or_whitespace_version_raises(self, version):
        with pytest.raises(ValueError, match="could not work out a name"):
            BuildManager.resolve_default_asset(
                "{{values.custom_seek_keys.version}}", self._namespaces(version)
            )

    @pytest.mark.parametrize(
        "template, version, expected",
        [
            # fasta/fasta_index declare `default_asset: "default"` by design.
            pytest.param("default", "", "default", id="literal-default-allowed"),
            pytest.param(
                "{{values.custom_seek_keys.version}}", "2.2.0", "2.2.0", id="real-version"
            ),
        ],
    )
    def test_nonempty_result_passes_through(self, template, version, expected):
        assert BuildManager.resolve_default_asset(template, self._namespaces(version)) == expected


PROBE = "refgenie.managers.build.check_output"


class TestCustomSeekKeyResolutionRunsWhereTheToolIs:
    """`BuildManager.resolve_custom_seek_keys`: where a version probe runs (unit tier).

    Asset names come from tool versions, and a tool version can only be learned
    by running the tool. A recipe that declares a `docker_image` is saying the
    tool lives in that image and not on this host -- 15 of the 29 stock recipes
    do exactly that. Probing the host anyway is why the web build form refused
    every version-named recipe: `bowtie2-build` is not installed on a server,
    correctly, so the name resolved empty and preflight blocked the build.

    The build path was always container-aware; this is the naming path being
    made to agree with it, through the same helper.
    """

    @staticmethod
    def _recipe(custom_seek_keys, docker_image=None, version="0.0.1"):
        """The three recipe attributes the resolver reads, and nothing else."""
        return SimpleNamespace(
            name="bowtie2_index",
            version=version,
            docker_image=docker_image,
            custom_seek_keys=custom_seek_keys,
        )

    def test_recipe_with_a_docker_image_probes_the_container(self, refgenie_minimal):
        recipe = self._recipe({"version": "bowtie2-build --version"}, "databio/refgenie")
        with patch(PROBE, return_value=b"2.3.0\n") as probe:
            assert refgenie_minimal.build.resolve_custom_seek_keys(recipe) == {"version": "2.3.0"}
        argv = probe.call_args.args[0]
        assert argv[:3] == ["docker", "run", "--rm"], argv
        assert "databio/refgenie" in argv
        assert argv[-1] == "bowtie2-build --version"

    def test_the_container_probe_is_one_shot_and_needs_no_tty(self, refgenie_minimal):
        """Not `docker run -itd` + `docker exec -it`: the detached container
        outlived every probe, and `exec -it` refuses whenever stdin is not a
        terminal -- a server thread, a snakemake run, CI -- printing to stderr
        and handing back an empty version."""
        recipe = self._recipe({"version": "bowtie2-build --version"}, "databio/refgenie")
        with patch(PROBE, return_value=b"2.3.0\n") as probe:
            refgenie_minimal.build.resolve_custom_seek_keys(recipe)
        assert probe.call_count == 1
        argv = probe.call_args.args[0]
        assert not any(flag in argv for flag in ("-it", "-itd", "exec"))
        assert "--rm" in argv

    def test_recipe_without_a_docker_image_still_probes_the_host(self, refgenie_minimal):
        """No image declared means the tool is expected on PATH, as before."""
        recipe = self._recipe({"version": "printf 1.2.3"})
        assert refgenie_minimal.build.resolve_custom_seek_keys(recipe) == {"version": "1.2.3"}

    def test_a_probe_that_prints_nothing_names_itself(self, refgenie_minimal):
        """The stock `bowtie2_index` case. `databio/refgenie` still ships
        bowtie2 2.3.0, which prints `bowtie2-build version 2.3.0`, while the
        recipe greps for `bowtie2-build-s version ` or `bowtie2-<digit>`. The
        grep matches nothing and still exits 0 through the pipe, so the probe
        succeeds with empty output. Surfaced only as an empty asset name, that
        told the user to check the commands run in this environment -- which
        they do. It has to name the command and the image instead.
        """
        recipe = self._recipe(
            {"version": "bowtie2-build --version | grep -oP x"}, "databio/refgenie"
        )
        with patch(PROBE, return_value=b"\n"):
            with pytest.raises(CustomSeekKeyError) as caught:
                refgenie_minimal.build.resolve_custom_seek_keys(recipe)
        message = str(caught.value)
        assert "names it after the version" in message
        assert "databio/refgenie" in message
        assert "printed nothing" in message
        assert "bowtie2-build --version" in message

    def test_a_recipe_with_no_custom_seek_keys_starts_nothing(self, refgenie_minimal):
        """Every recipe naming its asset the literal `default` is in this
        branch (abundant_sequences, bed12, blacklist, fasta...). Starting a
        container to name one of those would be pure cost."""
        recipe = self._recipe({}, "databio/refgenie")
        with patch(PROBE) as probe:
            assert refgenie_minimal.build.resolve_custom_seek_keys(recipe) == {}
        probe.assert_not_called()

    @pytest.mark.parametrize(
        "failure, expected",
        [
            pytest.param(
                FileNotFoundError(2, "No such file or directory: 'docker'"),
                "could not be started",
                id="docker-not-installed",
            ),
            pytest.param(
                CalledProcessError(125, "docker"),
                "exiting with status 125",
                id="daemon-down-or-no-image",
            ),
            pytest.param(
                TimeoutExpired("docker", CUSTOM_SEEK_KEY_TIMEOUT),
                "still running after",
                id="hung-probe-is-killed",
            ),
        ],
    )
    def test_unusable_docker_raises_a_named_error(self, refgenie_minimal, failure, expected):
        """Same as the build path: loud. A preflight catches it and renders it
        as a field problem, which is the honest answer -- far better than a
        wrong asset name baked into a published key."""
        recipe = self._recipe({"version": "bowtie2-build --version"}, "databio/refgenie")
        with patch(PROBE, side_effect=failure):
            with pytest.raises(CustomSeekKeyError) as raised:
                refgenie_minimal.build.resolve_custom_seek_keys(recipe)
        message = str(raised.value)
        assert expected in message, message
        assert "databio/refgenie" in message
        assert "version" in message

    def test_every_probe_is_bounded_by_a_timeout(self, refgenie_minimal):
        """`check_output` with no timeout on a `docker run` pins the request
        thread that called it for the life of the process."""
        for image in (None, "databio/refgenie"):
            with patch(PROBE, return_value=b"2.3.0\n") as probe:
                refgenie_minimal.build.resolve_custom_seek_keys(
                    self._recipe({"version": "bowtie2-build --version"}, image)
                )
            assert probe.call_args.kwargs["timeout"] == CUSTOM_SEEK_KEY_TIMEOUT

    def test_the_answer_is_cached_for_the_process(self, refgenie_minimal):
        """The build form preflights on a 400 ms debounce. A container per
        keystroke is not a price worth paying for a string that cannot change."""
        recipe = self._recipe({"version": "bowtie2-build --version"}, "databio/refgenie")
        with patch(PROBE, return_value=b"2.3.0\n") as probe:
            for _ in range(3):
                assert refgenie_minimal.build.resolve_custom_seek_keys(recipe) == {
                    "version": "2.3.0"
                }
        assert probe.call_count == 1

    def test_a_new_recipe_version_is_probed_again(self, refgenie_minimal):
        """The cache key is the recipe, not the command: a new recipe version
        may well name a new tool version."""
        with patch(PROBE, return_value=b"2.3.0\n") as probe:
            for version in ("0.0.1", "0.0.2"):
                refgenie_minimal.build.resolve_custom_seek_keys(
                    self._recipe({"version": "bowtie2-build --version"}, None, version)
                )
        assert probe.call_count == 2

    def test_a_failed_probe_is_not_cached(self, refgenie_minimal):
        """So `docker pull` and retry works, rather than needing a restart."""
        recipe = self._recipe({"version": "bowtie2-build --version"}, "databio/refgenie")
        with patch(PROBE, side_effect=FileNotFoundError(2, "docker")):
            with pytest.raises(CustomSeekKeyError):
                refgenie_minimal.build.resolve_custom_seek_keys(recipe)
        with patch(PROBE, return_value=b"2.3.0\n"):
            assert refgenie_minimal.build.resolve_custom_seek_keys(recipe) == {"version": "2.3.0"}


# --- Seek keys: builder, CLI parsing, persistence (unit tier) ----------------


class TestSeekKeyCLIParsing:
    """CLI --seek-key flag parsing via handle_add (unit tier)."""

    def _cmd(self, seek_keys):
        cmd = MagicMock()
        cmd.seek_keys = seek_keys
        cmd.asset_registry_paths = ["genome/group:asset"]
        cmd.genome_digest = None
        return cmd

    def _refgenie(self):
        r = MagicMock()
        r.parse_asset_registry_path.return_value = MagicMock(
            genome="genome", asset_group="group", asset="asset"
        )
        return r

    @pytest.mark.parametrize(
        "seek_keys,expected",
        [
            (["key=a=b"], {"key": "a=b"}),  # split on the first '=' only
            (["k1=v1", "k2=v2"], {"k1": "v1", "k2": "v2"}),
        ],
    )
    def test_valid_seek_keys(self, seek_keys, expected):
        r = self._refgenie()
        handle_add(self._cmd(seek_keys), r)
        assert r.asset.content.add.call_args[1]["custom_seek_keys"] == expected

    @pytest.mark.parametrize("seek_keys", [["no_equals_sign"], ["=value"]])
    def test_invalid_seek_key_format_raises(self, seek_keys):
        """A key with no '=' or an empty name raises ValueError."""
        with pytest.raises(ValueError, match="Invalid seek key format"):
            handle_add(self._cmd(seek_keys), MagicMock())

    def test_no_seek_keys_passes_none(self):
        """No --seek-key flag results in custom_seek_keys=None."""
        r = self._refgenie()
        handle_add(self._cmd(None), r)
        assert r.asset.content.add.call_args[1]["custom_seek_keys"] is None


class TestAssetClassSeekKeyPersistence:
    """Registering an asset-class YAML persists its seek keys; recipes are rejected (unit tier)."""

    def test_persists_prefix_seek_key(self, refgenie_minimal, fixtures_path):
        """A prefix-type seek key is persisted with its name, value, and type."""
        refgenie_minimal.asset_class.add(fixtures_path / "bowtie2_index_asset_class.yaml")

        retrieved = refgenie_minimal.asset_class.get("bowtie2_index", "0.0.1")
        assert len(retrieved.seek_keys) == 1
        sk = retrieved.seek_keys[0]
        assert sk.name == "bowtie2_index"
        assert sk.value == "{genome}"
        assert sk.type.value == "prefix"

    def test_persists_multiple_seek_keys(self, refgenie_minimal, fixtures_path):
        """An asset class with multiple seek keys persists all of them."""
        refgenie_minimal.asset_class.add(fixtures_path / "fasta_asset_class.yaml")

        retrieved = refgenie_minimal.asset_class.get("fasta", "0.1.0")
        assert {sk.name for sk in retrieved.seek_keys} == {"fasta", "fai", "chrom_sizes"}

    def test_recipe_yaml_rejected(self, refgenie_minimal, fixtures_path):
        """A recipe YAML passed to asset-class add is rejected, not silently accepted.

        Recipe YAMLs have output_asset_class/command_templates fields and no
        seek_keys; the code must detect this rather than create an empty class.
        """
        with pytest.raises(ValueError, match="recipe"):
            refgenie_minimal.asset_class.add(fixtures_path / "bowtie2_index_asset_recipe.yaml")


# --- merged from test_populate.py: populate/insert ---------------------------
# Populate replaces refgenie:// patterns in strings/files with resolved local
# paths. Insert (add) registers an externally-created asset directory into
# refgenie's database. Both classes are ``component``. The looper hook that
# builds on populate is tested in tests/integrations/test_looper.py.


class TestPopulate:
    """Test populate (refgenie:// path replacement)."""

    pytestmark = pytest.mark.component

    def test_populate_simple_string(self, refgenie_session):
        """Replace a single refgenie:// path in a string."""
        input_str = "samtools index refgenie://rCRSd/fasta:test"
        result = refgenie_session.populate(input_str)
        assert "refgenie://" not in result
        assert "/" in result

    def test_populate_preserves_non_refgenie_text(self, refgenie_session):
        """Non-refgenie text is preserved unchanged."""
        input_str = "echo hello world"
        result = refgenie_session.populate(input_str)
        assert result == input_str

    def test_populate_multiple_paths_in_string(self, refgenie_session):
        """Multiple refgenie:// patterns in one string are all replaced."""
        input_str = "cat refgenie://rCRSd/fasta:test | head refgenie://rCRSd/fasta:test"
        result = refgenie_session.populate(input_str)
        assert result.count("refgenie://") == 0

    def test_populate_nonexistent_genome_raises(self, refgenie_session):
        """Non-existent genome alias raises MissingAliasError."""
        input_str = "cat refgenie://nonexistent/fake:asset"
        with pytest.raises(MissingAliasError):
            refgenie_session.populate(input_str)

    def test_populate_missing_asset_is_left_as_written(self, refgenie_session):
        """A missing group on a known genome only warns; the text is unchanged."""
        input_str = "cat refgenie://rCRSd/nope"
        assert refgenie_session.populate(input_str) == input_str

    def test_populate_list_input(self, refgenie_session):
        """Populate works with list input."""
        input_list = [
            "line1 refgenie://rCRSd/fasta:test",
            "line2 no refgenie paths",
        ]
        result = refgenie_session.populate(input_list)
        assert isinstance(result, list)
        assert len(result) == 2
        assert "refgenie://" not in result[0]
        assert result[1] == "line2 no refgenie paths"

    def test_populate_dict_input(self, refgenie_session):
        """Populate works with dict input."""
        input_dict = {
            "key1": "refgenie://rCRSd/fasta:test",
            "key2": "no refgenie paths",
        }
        result = refgenie_session.populate(input_dict)
        assert isinstance(result, dict)
        assert "refgenie://" not in result["key1"]
        assert result["key2"] == "no refgenie paths"


class TestPopulateFile:
    """Test populate_file (file-based populate)."""

    pytestmark = pytest.mark.component

    def test_populate_file_outputs_to_stdout(self, refgenie_session, tmp_path):
        """populate_file reads a file and writes populated lines to stdout."""
        test_file = tmp_path / "test_config.txt"
        test_file.write_text("input: refgenie://rCRSd/fasta:test\noutput: /tmp/out\n")

        from refgenie.utils.io import populate_file

        captured = StringIO()
        with patch.object(sys, "stdout", captured):
            populate_file(file_path=test_file, pop_fun=refgenie_session.populate)

        output = captured.getvalue()
        assert "refgenie://" not in output
        assert "output: /tmp/out" in output


class TestInsert:
    """Test insert (adding external assets)."""

    pytestmark = pytest.mark.component

    def test_insert_new_asset(self, tmp_path, fixtures_path):
        """Insert an external asset directory into refgenie."""
        rg = Refgenie(database_engine=make_engine(), suppress_migrations=True)
        rg.database.init(genome_folder=tmp_path / "genomes")
        register_fasta(rg, fixtures_path)

        # Initialize genome first (need a genome to exist)
        rg.genome.initialize_genome(
            fasta_file_path=fixtures_path / "rCRSd.fa",
            alias_names=["rCRSd"],
            description="rCRSd genome",
        )
        rg.build.run(
            recipe_name="fasta",
            genome_alias="rCRSd",
            asset_group_name="fasta",
            asset_name="default",
        )

        # Create a fake asset directory inside genome_folder with the seek key files
        # the fasta asset class requires: {genome}.fa, {genome}.fa.fai, {genome}.chrom.sizes
        genome_digest = rg.alias.resolve("rCRSd")
        asset_dir = tmp_path / "genomes" / genome_digest / "custom_group" / "custom_asset"
        asset_dir.mkdir(parents=True)
        shutil.copy(fixtures_path / "rCRSd.fa", asset_dir / "rCRSd.fa")
        (asset_dir / "rCRSd.fa.fai").write_text("rCRSd\t16569\t7\t70\t71\n")
        (asset_dir / "rCRSd.chrom.sizes").write_text("rCRSd\t16569\n")

        result = rg.asset.content.add(
            asset_class_name="fasta",
            path=Path(genome_digest) / "custom_group" / "custom_asset",
            genome_digest=genome_digest,
            asset_group_name="custom_group",
            asset_name="custom_asset",
        )

        assert result is not None
        assert rg.asset.exists(
            genome_digest=genome_digest,
            asset_group_name="custom_group",
            asset_name="custom_asset",
        )

    def test_insert_nonexistent_path_raises(self, tmp_path, fixtures_path):
        """Insert with a non-existent path raises FileNotFoundError."""
        rg = Refgenie(database_engine=make_engine(), suppress_migrations=True)
        rg.database.init(genome_folder=tmp_path / "genomes")
        register_fasta(rg, fixtures_path)

        # Initialize genome first so the genome exists
        rg.genome.initialize_genome(
            fasta_file_path=fixtures_path / "rCRSd.fa",
            alias_names=["rCRSd"],
            description="rCRSd genome",
        )
        rg.build.run(
            recipe_name="fasta",
            genome_alias="rCRSd",
            asset_group_name="fasta",
            asset_name="default",
        )

        with pytest.raises(FileNotFoundError):
            rg.asset.content.add(
                asset_class_name="fasta",
                path=Path("nonexistent/path/to/asset"),
                genome_digest=rg.alias.resolve("rCRSd"),
                asset_group_name="fake_group",
                asset_name="fake_asset",
            )


# --- merged from test_seek_keys.py: non-path seek keys ------------------------
# Non-path (string/json) seek-key handling (seekr file mode is in
# test_servers.py). The class carries its own ``component`` mark: it builds real
# DB rows and genome folders. The unit half -- the builtin seek-key builder, CLI
# --seek-key parsing and asset-class seek-key persistence -- is above.


@pytest.fixture
def refgenie_with_test_asset_class(engine, tmp_path, fixtures_path):
    """Refgenie with the 'test' asset class (path + non-path seek keys) registered.

    Inits into tmp_path so these tests never write to the real
    ``config.genome_folder`` default.
    """
    r = Refgenie(database_engine=engine, suppress_migrations=True)
    r.database.init(
        genome_folder=tmp_path / "genomes",
        genome_stage_folder=tmp_path / "archives",
    )
    r.asset_class.add(fixtures_path / "test_asset_class.yaml")
    return r


# All non-path seek keys required by the "test" asset class.
ALL_NON_PATH_SEEK_KEYS = {
    "test_version": "1.2.3",
    "test_metadata": json.dumps({"key": "value"}),
}

RT_JSON = {"version": "1.0", "params": {"threads": 4}}


@pytest.fixture
def refgenie_test_class_genome(refgenie_with_test_asset_class, fixtures_path):
    """``refgenie_with_test_asset_class`` with the rCRSd genome initialized."""
    refgenie_with_test_asset_class.genome.initialize_genome(
        fasta_file_path=fixtures_path / "rCRSd.fa",
        alias_names=["rCRSd"],
        description="rCRSd genome",
    )
    return refgenie_with_test_asset_class


def _add_test_asset(r, asset_name, asset_dir_name, custom_seek_keys=None):
    """Materialize the three declared files and add a 'test'-class asset."""
    return add_asset_from_files(
        r,
        rel_dir=asset_dir_name,
        files=["{genome}.extension", "{genome}.extension1", "{genome}.extension2"],
        asset_class_name="test",
        asset_group_name="test",
        asset_name=asset_name,
        custom_seek_keys=custom_seek_keys or ALL_NON_PATH_SEEK_KEYS,
    )


class TestNonPathSeekKeys:
    """String/json/file seek keys through registration, seek, and servers.seek (component tier)."""

    pytestmark = pytest.mark.component

    def test_asset_class_registration_types(self, refgenie_with_test_asset_class):
        """Registration records file keys with a value and non-path keys with none."""
        ac = refgenie_with_test_asset_class.asset_class.get("test")
        assert len(ac.seek_keys) == 5
        by_name = {sk.name: sk for sk in ac.seek_keys}

        for name in ("extension", "extension1", "extension2"):
            assert by_name[name].value is not None
            assert by_name[name].type == SeekKeyType.file

        assert by_name["test_version"].value is None
        assert by_name["test_version"].type == SeekKeyType.string
        assert by_name["test_metadata"].value is None
        assert by_name["test_metadata"].type == SeekKeyType.json

    def test_recipe_custom_seek_keys_field(self, refgenie_with_test_asset_class, fixtures_path):
        """A recipe's custom_seek_keys field is populated from its YAML."""
        r = refgenie_with_test_asset_class
        register_fasta(r, fixtures_path)
        r.asset_class.add(fixtures_path / "bwa_index_asset_class.yaml")
        r.recipe.add(fixtures_path / "bwa_index_asset_recipe.yaml")

        recipe = r.recipe.get("bwa_index")
        assert recipe.custom_seek_keys is not None
        assert "bwa_version" in recipe.custom_seek_keys

    @pytest.mark.parametrize(
        "custom_seek_keys, seek_key, expected",
        [
            (
                {"test_version": "2.0.0", "test_metadata": json.dumps({"x": 1})},
                "test_version",
                "2.0.0",
            ),
            (
                {"test_version": "1.0", "test_metadata": json.dumps(RT_JSON)},
                "test_metadata",
                json.dumps(RT_JSON),
            ),
            (
                {"test_version": "1.0", "test_metadata": (json.dumps({"a": 1}), SeekKeyType.json)},
                "test_metadata",
                json.dumps({"a": 1}),
            ),
            (
                {**ALL_NON_PATH_SEEK_KEYS, "unlisted_key": "some_value"},
                "unlisted_key",
                "some_value",
            ),
        ],
        ids=["string-key", "json-key", "tuple-value-form", "undeclared-key"],
    )
    def test_seek_returns_non_path_value(
        self, refgenie_test_class_genome, custom_seek_keys, seek_key, expected
    ):
        r = refgenie_test_class_genome
        _add_test_asset(r, "nonpath_test", "test_asset_nonpath", custom_seek_keys)
        digest = r.alias.resolve("rCRSd")
        assert r.asset.seek(digest, "test", "nonpath_test", seek_key_name=seek_key) == expected

    def test_seek_returns_absolute_path_for_file_type(self, refgenie_test_class_genome):
        """seek() for a file-type seek key returns an absolute str path."""
        r = refgenie_test_class_genome
        _add_test_asset(r, "path_test", "test_asset_path")
        result = r.asset.seek(
            r.alias.resolve("rCRSd"), "test", "path_test", seek_key_name="extension"
        )
        assert isinstance(result, str)
        assert Path(result).is_absolute()

    def test_seek_keys_dict_excludes_non_path_keys(self, refgenie_test_class_genome):
        """Asset.seek_keys_dict contains only path-based keys, as Path values."""
        r = refgenie_test_class_genome
        _add_test_asset(r, "dict_test", "test_asset_dict")
        asset = r.asset.get(
            genome_digest=r.alias.resolve("rCRSd"), asset_group_name="test", asset_name="dict_test"
        )
        skd = asset.seek_keys_dict

        assert set(skd) >= {"extension", "extension1", "extension2"}
        assert "test_version" not in skd
        assert "test_metadata" not in skd
        assert all(isinstance(v, Path) for v in skd.values())

    def test_servers_seek_for_non_path_keys(self, refgenie_test_class_genome):
        """servers.seek() returns non-path values directly, with no server contact."""
        r = refgenie_test_class_genome
        original_json = {"tool": "bwa", "version": "0.7.17"}
        _add_test_asset(
            r,
            "remote_test",
            "test_asset_remote",
            {"test_version": "3.0.0", "test_metadata": json.dumps(original_json)},
        )
        assert (
            r.servers.seek(
                genome_digest=r.alias.resolve("rCRSd"),
                asset_group_name="test",
                asset_name="remote_test",
                seek_key="test_version",
            )
            == "3.0.0"
        )
        remote_json = r.servers.seek(
            genome_digest=r.alias.resolve("rCRSd"),
            asset_group_name="test",
            asset_name="remote_test",
            seek_key="test_metadata",
        )
        assert json.loads(remote_json) == original_json

    def test_metadata_only_asset_class(self, engine, tmp_path, fixtures_path):
        """An asset class with only non-path seek keys registers with all values None."""
        r = Refgenie(database_engine=engine, suppress_migrations=True)
        r.database.init(
            genome_folder=tmp_path / "genomes",
            genome_stage_folder=tmp_path / "archives",
        )
        r.asset_class.add(fixtures_path / "metadata_only_asset_class.yaml")

        ac = r.asset_class.get("metadata_only")
        assert len(ac.seek_keys) == 2
        for sk in ac.seek_keys:
            assert sk.value is None
            assert sk.type in (SeekKeyType.string, SeekKeyType.json)
