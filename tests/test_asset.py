"""
Tests for the asset manager CRUD/registry surface (refgenie.managers.asset).

Write-path correctness (add_from_path atomicity, write ordering) lives in
test_asset_content.py; removal lives in test_asset_removal.py. This file
covers registry-path parsing/population, the read-only query surface
(exists/seek/get/list/seek-keys/table), the AssetClass registration lifecycle,
and the unit half of seek-key handling (the built-in seek-keys builder, CLI
parsing, and asset-class seek-key persistence).

The component half of asset names lives in test_asset_content.py; the
component half of seek keys (non-path seek keys, seekr file mode) lives
below in this file.
"""

import json
import shutil
import sys
from io import StringIO
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from pydantic import ValidationError

from refgenie import Refgenie
from refgenie.cli.commands.curate import handle_add
from refgenie.db.tables import SeekKeyType
from refgenie.exceptions import (
    AssetClassExistsError,
    MissingAliasError,
    MissingAssetClassError,
    MissingAssetError,
)
from refgenie.managers.asset import AssetManager
from refgenie.populator import looper_refgenie_populate_local
from tests.helpers import (
    OMIT,
    add_asset_from_files,
    make_command_values,
    make_engine,
    mock_server_client,
    mocked_puller,
    register_fasta,
)


class TestAssetRegistryPath:
    """Parsing and population of refgenie:// registry paths (unit tier)."""

    @pytest.mark.parametrize(
        "asset_registry_path, result",
        [
            (
                "hg38/fasta:default",
                {"protocol": None, "genome": "hg38", "asset_group": "fasta",
                 "seek_key": None, "asset": "default"},
            ),
            (
                "hg38/fasta:custom",
                {"protocol": None, "genome": "hg38", "asset_group": "fasta",
                 "seek_key": None, "asset": "custom"},
            ),
            (
                "hg38/fasta.fasta:custom",
                {"protocol": None, "genome": "hg38", "asset_group": "fasta",
                 "seek_key": "fasta", "asset": "custom"},
            ),
            (
                "fasta",
                {"protocol": None, "genome": None, "asset_group": "fasta",
                 "seek_key": None, "asset": None},
            ),
            (
                "hg38/fasta",
                {"protocol": None, "genome": "hg38", "asset_group": "fasta",
                 "seek_key": None, "asset": None},
            ),
            (
                "test_proto://hg38/fasta",
                {"protocol": "test_proto", "genome": "hg38", "asset_group": "fasta",
                 "seek_key": None, "asset": None},
            ),
            (
                "refgenie://hg38/fasta.fasta:custom",
                {"protocol": "refgenie", "genome": "hg38", "asset_group": "fasta",
                 "seek_key": "fasta", "asset": "custom"},
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
        with patch.object(AssetManager, "_seek_by_components", return_value=Path("whatever")):
            assert (
                refgenie_minimal.populate_refgenie_registry_paths(asset_registry_path)
                == asset_registry_path
            )

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
        with patch.object(AssetManager, "_seek_by_components", return_value=str(Path(value))):
            assert refgenie_minimal.populate(asset_registry_path) == result

        with patch.object(AssetManager, "_seek_remote_by_components", return_value=value):
            assert refgenie_minimal.populater(asset_registry_path) == result


class TestAssetQuery:
    """Read-only asset query surface on a built catalog (unit tier)."""

    def test_exists(self, refgenie_session):
        """exists() is True for the built asset and False for each way to miss;
        get() returns the populated asset with its name and a non-empty digest."""
        r = refgenie_session
        assert r.asset.exists("fasta", "test", genome_name="rCRSd")
        assert not r.asset.exists("nonexistent", "test", genome_name="rCRSd")
        assert not r.asset.exists("fasta", "nonexistent", genome_name="rCRSd")
        assert not r.asset.exists(
            "fasta", "test", genome_digest="nonexistent_digest_123456789012345"
        )

        digest = r.alias.resolve("rCRSd")
        asset = r.asset.get(
            genome_digest=digest, asset_group_name="fasta", asset_name="test"
        )
        assert asset.name == "test"
        assert isinstance(asset.digest, str) and asset.digest

    def test_seek(self, refgenie_session):
        """seek returns a JSON-serializable str path to a real file; missing raises."""
        r = refgenie_session
        path = r.asset.seek("rCRSd", "fasta", "test")
        assert isinstance(path, str)
        json.dumps(path)  # seek values must be JSON-serializable
        assert Path(path).exists()
        # force_exists returns the same real path.
        assert r.asset.seek("rCRSd", "fasta", "test", force_exists=True) == path
        with pytest.raises(MissingAssetError):
            r.asset.seek("rCRSd", "nonexistent", "test", force_exists=True)

    def test_seek_accepts_genome_digest(self, refgenie_session):
        """seek(<digest>/fasta) resolves the same as seek(<alias>/fasta).

        The genome has a local alias, so a supplied digest falls back to that
        alias tree and yields an identical path.
        """
        r = refgenie_session
        digest = r.alias.resolve("rCRSd")
        assert digest != "rCRSd"
        by_alias = r.asset.seek("rCRSd", "fasta", "test")
        by_digest = r.asset.seek(digest, "fasta", "test")
        assert by_digest == by_alias

    def test_resolve_genome_digest_accepts_digest(self, refgenie_session):
        """The manager resolver accepts a known digest, and still rejects an
        unknown token as a MissingAliasError."""
        r = refgenie_session
        digest = r.alias.resolve("rCRSd")
        # A known digest passes through unchanged.
        assert r.asset._resolve_genome_digest(genome_digest=None, genome_name=digest) == digest
        # The alias still resolves to the same digest.
        assert r.asset._resolve_genome_digest(genome_digest=None, genome_name="rCRSd") == digest
        assert r.asset._resolve_genome_digests(genome_names=[digest]) == [digest]
        # An unknown token is neither alias nor digest.
        with pytest.raises(MissingAliasError):
            r.asset._resolve_genome_digest(genome_digest=None, genome_name="not-a-real-token")

    def test_list_assets_accepts_digest(self, refgenie_session):
        """list_assets (behind `list -g <digest>`) accepts a genome digest."""
        r = refgenie_session
        digest = r.alias.resolve("rCRSd")
        assets = list(r.asset.list_assets(genome_names=[digest]))
        assert assets
        assert "fasta" in {a.asset_group.name for a in assets}

    def test_list_assets_and_groups(self, refgenie_session):
        """list_assets returns the built asset(s), and the fasta group is among them."""
        assets = list(refgenie_session.asset.list_assets())
        assert assets
        assert "fasta" in {a.asset_group.name for a in assets}

    def test_list_seek_keys(self, refgenie_session):
        """list_seek_keys returns the fasta seek-key names; the default asset
        resolves to the same keys as naming it explicitly."""
        keys = refgenie_session.asset.list_seek_keys("rCRSd", "fasta")
        assert keys
        assert all(isinstance(k, str) for k in keys)
        assert "fasta" in keys
        assert keys == refgenie_session.asset.list_seek_keys(
            "rCRSd", "fasta", asset_name="test"
        )

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
        with pytest.raises(ValueError, match="empty asset name"):
            Refgenie.resolve_default_asset(
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
        assert Refgenie.resolve_default_asset(template, self._namespaces(version)) == expected


# --- Seek keys: builder, CLI parsing, persistence (unit tier) ----------------


class TestSeekKeyCLIParsing:
    """CLI --seek-key flag parsing via handle_add (unit tier)."""

    def _cmd(self, seek_keys):
        cmd = MagicMock()
        cmd.seek_keys = seek_keys
        cmd.asset_registry_paths = ["genome/group:asset"]
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
        assert r.add.call_args[1]["custom_seek_keys"] == expected

    @pytest.mark.parametrize("seek_keys", [["no_equals_sign"], ["=value"]])
    def test_invalid_seek_key_format_raises(self, seek_keys):
        """A key with no '=' or an empty name raises ValueError."""
        with pytest.raises(ValueError, match="Invalid seek key format"):
            handle_add(self._cmd(seek_keys), MagicMock())

    def test_no_seek_keys_passes_none(self):
        """No --seek-key flag results in custom_seek_keys=None."""
        r = self._refgenie()
        handle_add(self._cmd(None), r)
        assert r.add.call_args[1]["custom_seek_keys"] is None


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


# --- merged from test_populate.py: populate/insert + looper populator hook ---
# Populate replaces refgenie:// patterns in strings/files with resolved local
# paths. Insert (add) registers an externally-created asset directory into
# refgenie's database. The populator is the looper pre_submit hook that builds
# ``namespaces['refgenie'] = {genome: {asset_group: {seek_key: path}}}`` from a
# local refgenie1 SQLite database (drop-in for legacy
# ``refgenconf.looper_refgenie_populate``). Populate/insert classes are
# ``component``; the populator hook test is mock-light and stays ``unit``.


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
        rg.init(genome_folder=tmp_path / "genomes")
        register_fasta(rg, fixtures_path)

        # Initialize genome first (need a genome to exist)
        rg.genome.initialize_genome(
            fasta_file_path=fixtures_path / "rCRSd.fa",
            alias_names=["rCRSd"],
            description="rCRSd genome",
        )
        rg.build_asset(
            recipe_name="fasta",
            genome_name="rCRSd",
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

        result = rg.add(
            asset_class_name="fasta",
            path=Path(genome_digest) / "custom_group" / "custom_asset",
            genome_name="rCRSd",
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
        rg.init(genome_folder=tmp_path / "genomes")
        register_fasta(rg, fixtures_path)

        # Initialize genome first so the genome exists
        rg.genome.initialize_genome(
            fasta_file_path=fixtures_path / "rCRSd.fa",
            alias_names=["rCRSd"],
            description="rCRSd genome",
        )
        rg.build_asset(
            recipe_name="fasta",
            genome_name="rCRSd",
            asset_group_name="fasta",
            asset_name="default",
        )

        with pytest.raises(FileNotFoundError):
            rg.add(
                asset_class_name="fasta",
                path=Path("nonexistent/path/to/asset"),
                genome_name="rCRSd",
                asset_group_name="fake_group",
                asset_name="fake_asset",
            )


# --- the looper populator hook (unit) ---------------------------------------


def _empty_namespaces():
    return {
        "pipeline": {"var_templates": {}},
        "sample": {},
        "project": {},
    }


@pytest.fixture
def populated(refgenie_session, monkeypatch):
    """Populator output with populator.Refgenie() bound to the fixture session.

    Replaces the Refgenie() construction inside the hook with the pre-built,
    fixture-bound session instance, then runs the hook once on empty namespaces.
    """
    import refgenie.populator as populator

    class _Bound:
        def __init__(self, *args, **kwargs):
            # Discard args; use the fixture-built session instance
            self._inst = refgenie_session

        def __getattr__(self, item):
            return getattr(self._inst, item)

    monkeypatch.setattr(populator, "Refgenie", _Bound)
    return looper_refgenie_populate_local(_empty_namespaces())


def test_populator_seek_paths_match_seek_call(refgenie_session, populated):
    """The populator's emitted path equals what a direct seek() returns.

    Also asserts the emitted structure: a top-level 'refgenie' dict keyed by the
    genome alias, the asset-group key, and the fasta seek key -- all present with
    string leaves so jinja templating renders cleanly.
    """
    # Top-level shape: {'refgenie': {genome_alias: {asset_group: {seek_key: str}}}}
    assert isinstance(populated, dict)
    assert "refgenie" in populated
    assert isinstance(populated["refgenie"], dict)

    # Genome alias key (the session fixture builds a fasta asset on rCRSd)
    assert "rCRSd" in populated["refgenie"]
    rCRSd = populated["refgenie"]["rCRSd"]

    # Asset-group key and built-in fasta seek key
    assert "fasta" in rCRSd
    fasta_seek_keys = rCRSd["fasta"]
    assert "fasta" in fasta_seek_keys

    # Every leaf is a string (not Path) so jinja templating renders cleanly
    for sk_name, sk_val in fasta_seek_keys.items():
        assert isinstance(sk_val, str), f"seek_key {sk_name!r} value {sk_val!r} is not str"

    # The emitted path equals an authoritative seek() call.
    direct = str(refgenie_session.asset.seek("rCRSd", "fasta", "test", "fasta"))
    assert fasta_seek_keys["fasta"] == direct


# --- merged from test_seek_keys.py: non-path seek keys + seekr file mode ------
# Non-path (string/json) seek-key handling and seek_remote (seekr) file-mode URL
# construction. Both classes carry their own ``component`` mark: they build real
# DB rows and genome folders. The unit half -- the builtin seek-key builder, CLI
# --seek-key parsing and asset-class seek-key persistence -- is above.


@pytest.fixture
def refgenie_with_test_asset_class(engine, tmp_path, fixtures_path):
    """Refgenie with the 'test' asset class (path + non-path seek keys) registered.

    Inits into tmp_path so these tests never write to the real
    ``config.genome_folder`` default.
    """
    r = Refgenie(database_engine=engine, suppress_migrations=True)
    r.init(
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
    """String/json/file seek keys through registration, seek, and seek_remote (component tier)."""

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
            ({"test_version": "2.0.0", "test_metadata": json.dumps({"x": 1})},
             "test_version", "2.0.0"),
            ({"test_version": "1.0", "test_metadata": json.dumps(RT_JSON)},
             "test_metadata", json.dumps(RT_JSON)),
            ({"test_version": "1.0",
              "test_metadata": (json.dumps({"a": 1}), SeekKeyType.json)},
             "test_metadata", json.dumps({"a": 1})),
            ({**ALL_NON_PATH_SEEK_KEYS, "unlisted_key": "some_value"},
             "unlisted_key", "some_value"),
        ],
        ids=["string-key", "json-key", "tuple-value-form", "undeclared-key"],
    )
    def test_seek_returns_non_path_value(
        self, refgenie_test_class_genome, custom_seek_keys, seek_key, expected
    ):
        r = refgenie_test_class_genome
        _add_test_asset(r, "nonpath_test", "test_asset_nonpath", custom_seek_keys)
        assert r.asset.seek("rCRSd", "test", "nonpath_test", seek_key_name=seek_key) == expected

    def test_seek_returns_absolute_path_for_file_type(self, refgenie_test_class_genome):
        """seek() for a file-type seek key returns an absolute str path."""
        r = refgenie_test_class_genome
        _add_test_asset(r, "path_test", "test_asset_path")
        result = r.asset.seek("rCRSd", "test", "path_test", seek_key_name="extension")
        assert isinstance(result, str)
        assert Path(result).is_absolute()

    def test_seek_keys_dict_excludes_non_path_keys(self, refgenie_test_class_genome):
        """Asset.seek_keys_dict contains only path-based keys, as Path values."""
        r = refgenie_test_class_genome
        _add_test_asset(r, "dict_test", "test_asset_dict")
        asset = r.asset.get(genome_name="rCRSd", asset_group_name="test", asset_name="dict_test")
        skd = asset.seek_keys_dict

        assert set(skd) >= {"extension", "extension1", "extension2"}
        assert "test_version" not in skd
        assert "test_metadata" not in skd
        assert all(isinstance(v, Path) for v in skd.values())

    def test_seek_remote_for_non_path_keys(self, refgenie_test_class_genome):
        """seek_remote() returns non-path values directly, with no server contact."""
        r = refgenie_test_class_genome
        original_json = {"tool": "bwa", "version": "0.7.17"}
        _add_test_asset(
            r, "remote_test", "test_asset_remote",
            {"test_version": "3.0.0", "test_metadata": json.dumps(original_json)},
        )
        assert r.asset.seek_remote(
            genome_name="rCRSd", asset_group_name="test",
            asset_name="remote_test", seek_key="test_version",
        ) == "3.0.0"
        remote_json = r.asset.seek_remote(
            genome_name="rCRSd", asset_group_name="test",
            asset_name="remote_test", seek_key="test_metadata",
        )
        assert json.loads(remote_json) == original_json

    def test_metadata_only_asset_class(self, engine, tmp_path, fixtures_path):
        """An asset class with only non-path seek keys registers with all values None."""
        r = Refgenie(database_engine=engine, suppress_migrations=True)
        r.init(
            genome_folder=tmp_path / "genomes",
            genome_stage_folder=tmp_path / "archives",
        )
        r.asset_class.add(fixtures_path / "metadata_only_asset_class.yaml")

        ac = r.asset_class.get("metadata_only")
        assert len(ac.seek_keys) == 2
        for sk in ac.seek_keys:
            assert sk.value is None
            assert sk.type in (SeekKeyType.string, SeekKeyType.json)


# --- seekr (seek_remote) file-mode URL construction ---------------------------


def _seekr_client(serving_modes, **kwargs):
    """A mock ServerClient for seekr testing (no genome_digest on the group)."""
    return mock_server_client(
        asset_group_name="fasta",
        genome_digest=OMIT,
        asset_name="test",
        asset_digest="asset_digest_001",
        serving_modes=serving_modes,
        seek_keys=[
            {"name": "fasta", "value": "rCRSd.fa", "type": "file"},
            {"name": "fai", "value": "rCRSd.fa.fai", "type": "file"},
        ],
        **kwargs,
    )


class TestSeekrFileMode:
    """seekr (seek remote) with file-mode assets, against the session-scoped
    built FASTA and a mocked server client (component tier)."""

    pytestmark = pytest.mark.component

    @pytest.mark.parametrize(
        "seek_key, expect_in_result",
        [
            (None, None),  # default seek key -> just the file endpoint
            ("fai", ".fai"),  # named seek key -> URL references that seek key's file
        ],
    )
    def test_seekr_returns_file_url_for_file_mode_asset(
        self, refgenie_session, seek_key, expect_in_result
    ):
        """seekr returns a file-endpoint URL for file-mode assets, using the requested
        seek key's file path when one is given."""
        r = refgenie_session
        mock_client = _seekr_client(serving_modes=["file"])

        kwargs = {} if seek_key is None else {"seek_key": seek_key}
        with mocked_puller(
            r,
            mock_client,
            mock_genome=False,
            mock_download_modes=False,
            mock_asset_writes=False,
        ):
            result = r.asset.seek_remote(
                genome_name="rCRSd",
                asset_group_name="fasta",
                asset_name="test",
                **kwargs,
            )

        assert result.startswith("http://test.example.com/v4/assets/asset_digest_001/files/")
        if expect_in_result is not None:
            assert expect_in_result in result

    def test_seekr_resolves_alias_via_server_without_local_genome(self, refgenie_session):
        """seekr answers for a genome the client does NOT have locally: the alias
        is resolved read-only against the server, and no local genome/alias row
        is created."""
        r = refgenie_session
        remote_digest = "0" * 32
        assert not r.alias.exists("notlocal")
        assert not r.genome.exists(remote_digest)

        mock_client = _seekr_client(serving_modes=["file"], alias_digest=remote_digest)
        with mocked_puller(
            r,
            mock_client,
            mock_genome=False,
            mock_download_modes=False,
            mock_asset_writes=False,
        ):
            result = r.asset.seek_remote(
                genome_name="notlocal",
                asset_group_name="fasta",
                asset_name="test",
            )

        assert result.startswith("http://test.example.com/v4/assets/asset_digest_001/files/")
        # Read-only: nothing was written locally.
        assert not r.alias.exists("notlocal")
        assert not r.genome.exists(remote_digest)

    def test_seekr_resolves_default_asset_from_server(self, refgenie_session):
        """With the asset name omitted, seekr picks the server's is_default asset
        for a non-local genome."""
        r = refgenie_session
        remote_digest = "1" * 32
        mock_client = _seekr_client(
            serving_modes=["file"], alias_digest=remote_digest, is_default=True
        )
        with mocked_puller(
            r,
            mock_client,
            mock_genome=False,
            mock_download_modes=False,
            mock_asset_writes=False,
        ):
            result = r.asset.seek_remote(
                genome_name="notlocal",
                asset_group_name="fasta",
            )
        assert result.startswith("http://test.example.com/v4/assets/asset_digest_001/files/")

    def test_seekr_errors_for_archive_only_asset(self, refgenie_session):
        """seekr raises ValueError for archive-only assets."""
        r = refgenie_session
        mock_client = _seekr_client(serving_modes=["archive"])

        with mocked_puller(
            r,
            mock_client,
            mock_genome=False,
            mock_download_modes=False,
            mock_asset_writes=False,
        ):
            with pytest.raises(ValueError, match="does not support file-level access"):
                r.asset.seek_remote(
                    genome_name="rCRSd",
                    asset_group_name="fasta",
                    asset_name="test",
                )
