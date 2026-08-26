"""CLI tests for refgenie: argument parsing, --help output, and dispatch.

All tests go through the pydantic-settings entry point
(``refgenie.cli.main.main``) or its parsing internals. The file is organized
in three groups:

1. **argv -> pydantic parsing and model validation** -- pure parsing tests
   that construct pydantic-settings models from argv without touching the
   database or filesystem.
2. **--help output** -- top-level command grouping, placeholder marking, and
   per-command help rendering.
3. **Dispatch and exit codes** -- commands run against a Refgenie instance
   backed by an in-memory SQLite engine patched onto
   ``Refgenie.get_default_database_engine``, asserting on exit codes, output,
   and dispatch (i.e. that the CLI wires up to the library correctly).

Classes that build a real FASTA asset on disk are marked ``component``
(deselected from the bare ``pytest`` loop); everything else is the default
``unit`` tier.
"""

import io
import sys
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from pydantic import AliasChoices, Field, ValidationError
from pydantic_settings import get_subcommand

from refgenie import Refgenie
from refgenie.cli import main as cli_main
from refgenie.cli.main import _preprocess_argv, main
from refgenie.cli.commands.alias import AliasGetNestedModel
from refgenie.cli.commands.data_channel import DataChannelSyncNestedModel
from refgenie.cli.commands.pull import MirrorModel, PullModel
from tests.helpers import make_engine, register_fasta


# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------


def parse(argv):
    """Parse argv through the pydantic-settings TopLevelParser."""
    from refgenie.cli.commands.framework import CleanHelpCliSource
    from refgenie.cli.parser import TopLevelParser

    source = CleanHelpCliSource(TopLevelParser, cli_parse_args=argv)
    return TopLevelParser(_cli_settings_source=source)


def subcmd(argv):
    """Parse argv and return the active leaf subcommand model."""
    return get_subcommand(parse(argv), is_required=True)


@pytest.fixture
def cli_engine():
    """In-memory SQLite engine for CLI tests."""
    return make_engine()


@pytest.fixture
def cli_rg(cli_engine, tmp_path):
    """Initialized Refgenie instance sharing the CLI engine."""
    rg = Refgenie(database_engine=cli_engine, suppress_migrations=True)
    rg.init(genome_folder=tmp_path / "genomes")
    return rg


@pytest.fixture
def cli_rg_with_recipe(tmp_path):
    """Initialized Refgenie with its own engine (used by build --requirements)."""
    rg = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    rg.init(genome_folder=tmp_path / "genomes")
    return rg


@pytest.fixture
def cli(cli_engine):
    """CLI invoker: invoke(*args) -> exit_code (use capsys to capture output)."""

    def invoke(*args):
        with patch.object(
            Refgenie, "get_default_database_engine", return_value=cli_engine
        ):
            import pydantic

            try:
                main(test_args=list(args))
                return 0
            except SystemExit as e:
                return e.code or 0
            except pydantic.ValidationError:
                return 2

    return invoke


@pytest.fixture
def run_main(cli_engine):
    """Run main() with the patched engine, letting SystemExit propagate."""

    def _run(*args):
        with patch.object(
            Refgenie, "get_default_database_engine", return_value=cli_engine
        ):
            main(test_args=list(args))

    return _run


# ===========================================================================
# argv -> pydantic parsing and model validation (no database, no filesystem)
# ===========================================================================


# ---------------------------------------------------------------------------
# Top-level and nested subcommand argument parsing
# ---------------------------------------------------------------------------


class TestPhase2Parsing:
    def test_subcommands_parse(self):
        """Each top-level subcommand routes to its model with the right fields."""
        assert parse(["seek", "hg38/fasta"]).seek.asset_registry_paths == ["hg38/fasta"]

        cmp = parse(["compare", "hg38", "mm10"]).compare
        assert cmp.genome1 == "hg38" and cmp.genome2 == "mm10"

        pull = parse(["pull", "--all", "-g", "hg38"]).pull
        assert getattr(pull, "all") is True and pull.genome == "hg38"
        # asset_registry_paths is non-positional with default None for `pull`
        # (pydantic-settings disallows positional args with defaults).
        assert pull.asset_registry_paths in (None, [])

        rm = parse(["remove", "hg38/fasta", "-f"]).remove
        assert rm.force is True and rm.asset_registry_paths == ["hg38/fasta"]

        assert parse(["listr", "-s", "http://example.com"]).listr.genome_server == [
            "http://example.com"
        ]
        assert parse(["mirror", "-f"]).mirror.force is True
        assert parse(["id", "hg38"]).id.asset_registry_paths == ["hg38"]
        assert (
            parse(["rename", "hg38/fasta", "-n", "fasta_v2"]).rename.new_asset_name == "fasta_v2"
        )

        gs = parse(["getseq", "-g", "hg38", "-l", "chr1:1-100"]).getseq
        assert gs.genome == "hg38" and gs.locus == "chr1:1-100"

        add = parse(["add", "hg38/fasta:custom", "-p", "/tmp/path", "-c", "fasta"]).add
        assert add.asset_registry_paths == ["hg38/fasta:custom"]
        assert add.path == "/tmp/path" and add.asset_class == "fasta"

        assert parse(["populate", "-f", "/tmp/file.txt"]).populate.file == "/tmp/file.txt"

        seekr = parse(["seekr", "hg38/fasta", "-s", "http://example.com"]).seekr
        assert seekr.asset_registry_paths == ["hg38/fasta"]
        assert seekr.genome_server == ["http://example.com"]


class TestPhase4NestedParsing:
    def test_nested_subcommands_parse(self):
        """Each nested subcommand routes to its model with the right fields."""
        ag = parse(["alias", "get", "-a", "hg38"]).alias.get
        assert ag.aliases == ["hg38"] and ag.genome_digests is None

        aset = parse(["alias", "set", "-a", "my_alias", "-d", "abc123"]).alias.set
        assert aset.aliases == ["my_alias"] and aset.digest == "abc123"

        assert parse(["config", "get"]).config.get is not None
        assert parse(["recipe", "list"]).recipe.list is not None
        assert parse(["recipe", "show", "fasta"]).recipe.show.recipe_name == "fasta"
        assert parse(["recipe", "add", "--source", "/path/to/recipe"]).recipe.add.source == (
            "/path/to/recipe"
        )
        assert parse(["asset-class", "list"]).asset_class.list is not None
        assert (
            parse(["asset-class", "show", "fasta"]).asset_class.show.asset_class_name == "fasta"
        )
        assert parse(["stage", "list"]).stage.list is not None

        dc_add = parse(
            ["data-channel", "add", "test_ch", "http", "http://x.com/index.yaml"]
        ).data_channel.add
        assert dc_add.name == "test_ch" and dc_add.type == "http"
        assert dc_add.index_address == "http://x.com/index.yaml"

        dc_sync = parse(["data-channel", "sync", "my_channel", "--exists-ok"]).data_channel.sync
        assert dc_sync.name == "my_channel"
        assert dc_sync.exists_ok is True and dc_sync.exists_overwrite is False

        snake = parse(["generate", "snakefile", "-o", "/tmp/Snakefile"]).generate.snakefile
        assert snake.output_path == Path("/tmp/Snakefile")

        r_add = parse(
            ["remote", "add", "--type", "http", "--prefix", "http://x.com", "--description", "test"]
        ).remote.add
        assert r_add.type == "http" and r_add.prefix == "http://x.com"
        assert r_add.description == "test"

        assert parse(["genome", "list"]).genome.list is not None

        g_rm = parse(["genome", "remove", "--genome", "rCRSd", "-f"]).genome.remove
        assert g_rm.genome == ["rCRSd"] and g_rm.force is True


# ---------------------------------------------------------------------------
# build command flag parsing
# ---------------------------------------------------------------------------


class TestBuildFlagParsing:
    def test_docker_set_true(self):
        assert subcmd(["build", "rCRSd/fasta", "--docker"]).docker is True

    def test_files_flat_list(self):
        args = subcmd(["build", "rCRSd/fasta", "--files", "fasta=/path/to/file.fa"])
        assert args.files == ["fasta=/path/to/file.fa"]

    def test_files_repeated_flag(self):
        args = subcmd(
            ["build", "rCRSd/fasta", "--files", "fasta=/a.fa", "--files", "gtf=/b.gtf"]
        )
        assert args.files == ["fasta=/a.fa", "gtf=/b.gtf"]


# ---------------------------------------------------------------------------
# Mutually exclusive groups -- validated at both the model and CLI-parse level
# ---------------------------------------------------------------------------


class TestPullModelValidation:
    """PullModel mutually exclusive group validators (direct construction)."""

    def test_skip_large_and_pull_large_raises(self):
        with pytest.raises(ValidationError, match="mutually exclusive"):
            PullModel(skip_large=True, pull_large=True)

    def test_skip_large_alone_ok(self):
        m = PullModel(skip_large=True)
        assert m.skip_large is True
        assert m.pull_large is False

    def test_pull_large_alone_ok(self):
        m = PullModel(pull_large=True)
        assert m.pull_large is True
        assert m.skip_large is False

    def test_batch_sets_pull_large(self):
        assert PullModel(batch=True).pull_large is True

    def test_batch_with_skip_large_raises(self):
        with pytest.raises(ValidationError, match="--batch.*conflicts.*--skip-large"):
            PullModel(batch=True, skip_large=True)

    def test_resolve_force_large_skip_large(self):
        assert PullModel(skip_large=True).resolve_force_large() is False

    def test_resolve_force_large_pull_large(self):
        assert PullModel(pull_large=True).resolve_force_large() is True

    def test_resolve_force_large_default(self):
        assert PullModel().resolve_force_large() is None


class TestMirrorModelValidation:
    """MirrorModel shares PullModel's validators via the mixin. One smoke test
    proves the mixin is wired up; the full matrix lives in
    TestPullModelValidation."""

    def test_skip_large_and_pull_large_raises(self):
        with pytest.raises(ValidationError, match="mutually exclusive"):
            MirrorModel(skip_large=True, pull_large=True)


class TestAliasGetNestedModelValidation:
    def test_aliases_and_genome_digests_raises(self):
        with pytest.raises(ValidationError, match="mutually exclusive"):
            AliasGetNestedModel(aliases=["hg38"], genome_digests=["abc123"])

    def test_aliases_alone_ok(self):
        m = AliasGetNestedModel(aliases=["hg38", "hg19"])
        assert m.aliases == ["hg38", "hg19"]
        assert m.genome_digests is None

    def test_genome_digests_alone_ok(self):
        m = AliasGetNestedModel(genome_digests=["abc123"])
        assert m.genome_digests == ["abc123"]
        assert m.aliases is None


class TestDataChannelSyncNestedModelValidation:
    def test_exists_ok_and_exists_overwrite_raises(self):
        with pytest.raises(ValidationError, match="mutually exclusive"):
            DataChannelSyncNestedModel(name="test_channel", exists_ok=True, exists_overwrite=True)

    def test_exists_ok_alone_ok(self):
        m = DataChannelSyncNestedModel(name="test_channel", exists_ok=True)
        assert m.exists_ok is True
        assert m.exists_overwrite is False

    def test_exists_overwrite_alone_ok(self):
        m = DataChannelSyncNestedModel(name="test_channel", exists_overwrite=True)
        assert m.exists_overwrite is True
        assert m.exists_ok is False


class TestMutualExclusionThroughCliParse:
    """The validators must also fire through the CLI parse path, not only
    direct model construction."""

    def test_pull_skip_large_pull_large(self):
        with pytest.raises(ValidationError, match="mutually exclusive"):
            parse(["pull", "--skip-large", "--pull-large"])

    def test_pull_batch_skip_large(self):
        with pytest.raises(ValidationError, match="conflicts with --skip-large"):
            parse(["pull", "--batch", "--skip-large"])

    def test_mirror_skip_large_pull_large(self):
        with pytest.raises(ValidationError, match="mutually exclusive"):
            parse(["mirror", "--skip-large", "--pull-large"])

    def test_data_channel_sync_exists_flags(self):
        with pytest.raises(ValidationError, match="mutually exclusive"):
            parse(
                ["data-channel", "sync", "test_ch", "--exists-ok", "--exists-overwrite"]
            )

    def test_alias_get_aliases_and_digests(self):
        with pytest.raises(ValidationError, match="mutually exclusive"):
            parse(["alias", "get", "-a", "x", "-g", "y"])

    def test_pull_size_cutoff(self):
        assert parse(["pull", "--size-cutoff", "5.5"]).pull.size_cutoff == 5.5

    def test_pull_batch_mode(self):
        args = parse(["pull", "--batch"])
        assert args.pull.pull_large is True
        assert args.pull.batch is True


class TestPullPositionalArgParsing:
    """Positional asset_registry_paths in pull command (via argv preprocessing)."""

    def _parse_pull(self, args_list):
        from refgenie.cli.main import _make_cli_source
        from refgenie.cli.parser import TopLevelParser

        argv = _preprocess_argv(["pull"] + args_list)
        cli_source = _make_cli_source(TopLevelParser, argv)
        return get_subcommand(TopLevelParser(_cli_settings_source=cli_source), is_required=True)

    def test_single_positional_path(self):
        assert self._parse_pull(["GRCm39/fasta"]).asset_registry_paths == ["GRCm39/fasta"]

    def test_multiple_positional_paths(self):
        assert self._parse_pull(["GRCm39/fasta", "GRCh38/fasta"]).asset_registry_paths == [
            "GRCm39/fasta",
            "GRCh38/fasta",
        ]

    def test_no_positional_args_with_flags(self):
        subcommand = self._parse_pull(["--all", "--genome", "GRCm39"])
        assert subcommand.asset_registry_paths in (None, [])
        assert getattr(subcommand, "all") is True
        assert subcommand.genome == "GRCm39"

    def test_positional_with_flags(self):
        subcommand = self._parse_pull(["GRCm39/fasta", "--force"])
        assert subcommand.asset_registry_paths == ["GRCm39/fasta"]
        assert subcommand.force is True


def _parse_pydantic_args(args_list):
    """Parse CLI args with pydantic-settings; return the subcommand model.

    Applies the same argv preprocessing as main() so that bare positional paths
    after ``pull`` are routed to ``--asset-registry-paths``.
    """
    from pydantic_settings import get_subcommand

    from refgenie.cli.main import _make_cli_source, _preprocess_argv
    from refgenie.cli.parser import TopLevelParser

    cli_source = _make_cli_source(TopLevelParser, _preprocess_argv(args_list))
    parsed = TopLevelParser(_cli_settings_source=cli_source)
    return get_subcommand(parsed, is_required=True)


class TestBulkPullArgParsing:
    """CLI argument parsing for bulk pull modes."""

    @pytest.mark.parametrize(
        "argv, model, expected",
        [
            (["pull", "-g", "hg38,mm10", "--all"], PullModel, {"genome": "hg38,mm10", "all": True}),
            (["pull", "--all-genomes", "--asset", "fasta"], PullModel,
             {"all_genomes": True, "asset": "fasta"}),
            (["pull", "-g", "hg38", "--init"], PullModel, {"genome": "hg38", "init": True}),
            (["pull", "-g", "hg38", "--all", "--force"], PullModel, {"force": True}),
            (["pull", "hg38/fasta"], PullModel,
             {"asset_registry_paths": ["hg38/fasta"], "all": False, "init": False}),
            (["mirror"], MirrorModel, {}),
            (["mirror", "--force"], MirrorModel, {"force": True}),
            (["mirror", "--size-cutoff", "20"], MirrorModel, {"size_cutoff": 20.0}),
        ],
        ids=[
            "pull-all-comma-genomes", "pull-all-genomes-asset", "pull-init",
            "pull-all-force", "pull-positional-path", "mirror", "mirror-force",
            "mirror-size-cutoff",
        ],
    )
    def test_bulk_pull_arg_parsing(self, argv, model, expected):
        subcmd = _parse_pydantic_args(argv)
        assert isinstance(subcmd, model)
        for name, value in expected.items():
            assert getattr(subcmd, name) == value, f"{name} mismatch for {' '.join(argv)}"


# ---------------------------------------------------------------------------
# argv preprocessing in refgenie.cli.main
# ---------------------------------------------------------------------------


class TestDerivedValueConsumingFlags:
    """`refgenie pull <paths>` is rescued by rewriting argv, because
    pydantic-settings rejects positional args that have defaults. The rewrite
    must learn which options consume a following value from the model, not a
    hand-maintained list."""

    def test_new_value_taking_option_is_handled(self, monkeypatch):
        """A value-taking option added to PullModel must not swallow its value
        as a registry path."""
        from refgenie.cli.commands import pull as commands_mod

        class PullWithTag(commands_mod.PullModel):
            tag: str | None = Field(None, validation_alias=AliasChoices("tag"))

        monkeypatch.setattr(commands_mod, "PullModel", PullWithTag)

        out = _preprocess_argv(["pull", "--tag", "default", "hg38/fasta"])
        idx = out.index("--asset-registry-paths")
        assert out[idx + 1] == "hg38/fasta"


class TestDataChannelAliasRewrite:
    def test_subcommand_position_is_rewritten(self):
        assert _preprocess_argv(["data_channel", "list"])[0] == "data-channel"

    def test_channel_named_data_channel_is_left_alone(self):
        out = _preprocess_argv(["data-channel", "show", "data_channel"])
        assert out == ["data-channel", "show", "data_channel"]


class TestAssetListAliasesTopLevelList:
    """`refgenie asset list` forwards to handle_list, so it must deserialize
    `-g a,b` the same way `refgenie list` does."""

    def test_comma_separated_genomes_match(self):
        top = parse(["list", "-g", "a,b"]).list
        nested = parse(["asset", "list", "-g", "a,b"]).asset.list
        assert nested.genome == top.genome == ["a", "b"]

    def test_genome_field_deserializes_like_list_model(self):
        """The two models must accept the same inputs, not just the same argv.

        Through argparse both already worked; direct construction did not --
        AssetListNestedModel.genome was list[str] while ListModel.genome is
        CliList, so only the latter accepted a comma-separated string.
        """
        from refgenie.cli.commands.asset import AssetListNestedModel
        from refgenie.cli.commands.listing import ListModel

        assert AssetListNestedModel(genome="a,b").genome == ListModel(genome="a,b").genome

    def test_asset_list_subcommand_routes_to_nested_model(self):
        """``refgenie asset list`` parses to AssetListNestedModel."""
        from refgenie.cli.main import _make_cli_source
        from refgenie.cli.commands.asset import AssetListNestedModel
        from refgenie.cli.parser import TopLevelParser

        src = _make_cli_source(TopLevelParser, ["asset", "list"])
        parsed = TopLevelParser(_cli_settings_source=src)
        asset_group = get_subcommand(parsed, is_required=True)
        leaf = get_subcommand(asset_group, is_required=True)
        assert isinstance(leaf, AssetListNestedModel)


class TestVersionFastPath:
    def test_version_as_first_arg_prints_version(self, capsys):
        with pytest.raises(SystemExit) as exc:
            cli_main.main(test_args=["--version"])
        assert exc.value.code == 0
        assert "refgenie" in capsys.readouterr().out

    def test_version_short_flag(self, capsys):
        with pytest.raises(SystemExit) as exc:
            cli_main.main(test_args=["-V"])
        out = capsys.readouterr().out
        assert exc.value.code == 0
        assert "refgenie" in out and any(ch.isdigit() for ch in out)

    def test_version_elsewhere_is_not_a_fast_path(self, capsys):
        """`refgenie seek --version` must not be treated as `refgenie --version`."""
        with pytest.raises(SystemExit):
            cli_main.main(test_args=["seek", "--version"])
        out = capsys.readouterr().out
        assert not out.startswith("refgenie ")


# ===========================================================================
# --help behavior
# ===========================================================================


class TestTopLevelHelp:
    """`refgenie --help` groups commands under headers from COMMAND_GROUPS."""

    @pytest.fixture(scope="class")
    def help_output(self):
        buf = io.StringIO()
        with redirect_stdout(buf):
            with pytest.raises(SystemExit):
                main(["--help"])
        return buf.getvalue()

    def test_all_commands_present(self, help_output):
        from refgenie.cli.messages import COMMAND_GROUPS

        for commands in COMMAND_GROUPS.values():
            for cmd in commands:
                assert cmd in help_output, f"Command '{cmd}' missing from --help"

    def test_group_headers_before_commands(self, help_output):
        """Each group header appears before all its member commands in the
        grouped section (searched after the usage/options block)."""
        import re

        from refgenie.cli.messages import COMMAND_GROUPS

        first_group = next(iter(COMMAND_GROUPS))
        grouped_section = help_output[help_output.index(first_group) :]

        for group_name, commands in COMMAND_GROUPS.items():
            header_pos = grouped_section.index(group_name)
            for cmd in commands:
                match = re.search(rf"^\s+{re.escape(cmd)}\s", grouped_section, re.MULTILINE)
                assert match is not None, f"Command '{cmd}' not found in grouped section"
                assert header_pos < match.start(), (
                    f"Group header '{group_name}' not before command '{cmd}'"
                )

    def test_no_flat_subcommands_block(self, help_output):
        assert "subcommands:" not in help_output


class TestPlaceholderHelp:
    """Commands that always fail must say so in --help, not only a docstring."""

    @pytest.mark.parametrize(
        "argv, subcommands",
        [
            (["config", "--help"], ["set"]),
        ],
    )
    def test_placeholder_marked(self, argv, subcommands, capsys):
        with pytest.raises(SystemExit):
            main(test_args=argv)
        out = capsys.readouterr().out
        for name in subcommands:
            line = next((ln for ln in out.splitlines() if ln.strip().startswith(name)), None)
            assert line is not None, f"{name} missing from help:\n{out}"
            assert "[not implemented]" in line, f"{name} not marked: {line!r}"


class TestHelpTextRendering:
    def test_no_build_help_uses_short_form(self, capsys):
        """`refgenie genome init --help` shows `--no-build`, not the dotted
        `genome.init.no-build`."""
        with pytest.raises(SystemExit):
            cli_main.main(test_args=["genome", "init", "--help"])
        captured = capsys.readouterr()
        text = captured.out + captured.err
        assert "--no-build" in text
        assert "genome.init.no-build" not in text


class TestSubcommandHelp:
    """Every command's --help exits 0 and produces output.

    Top-level ``--help`` and the ``config``/``recipe``/``asset-class`` group
    help are covered with stronger assertions by TestTopLevelHelp and
    TestPlaceholderHelp, so they are not repeated here.
    """

    @pytest.mark.parametrize(
        "cmd",
        [
            ["init", "--help"],
            ["purge", "--help"],
            ["list", "--help"],
            ["subscribe", "--help"],
            ["unsubscribe", "--help"],
            ["seek", "--help"],
            ["seekr", "--help"],
            ["remove", "--help"],
            ["rename", "--help"],
            ["id", "--help"],
            ["add", "--help"],
            ["getseq", "--help"],
            ["pull", "--help"],
            ["listr", "--help"],
            ["compare", "--help"],
            ["populate", "--help"],
            ["populater", "--help"],
            ["mirror", "--help"],
            ["build", "--help"],
            ["alias", "--help"],
            ["alias", "get", "--help"],
            ["alias", "set", "--help"],
            ["alias", "remove", "--help"],
            ["config", "get", "--help"],
            ["config", "set", "--help"],
            ["recipe", "list", "--help"],
            ["recipe", "show", "--help"],
            ["recipe", "add", "--help"],
            ["recipe", "remove", "--help"],
            ["recipe", "requirements", "--help"],
            ["asset-class", "list", "--help"],
            ["asset-class", "show", "--help"],
            ["asset-class", "add", "--help"],
            ["asset-class", "remove", "--help"],
            ["stage", "--help"],
            ["stage", "add", "--help"],
            ["stage", "remove", "--help"],
            ["stage", "list", "--help"],
            ["data-channel", "--help"],
            ["data-channel", "list", "--help"],
            ["data-channel", "add", "--help"],
            ["data-channel", "remove", "--help"],
            ["data-channel", "show", "--help"],
            ["data-channel", "validate", "--help"],
            ["data-channel", "sync", "--help"],
            ["generate", "--help"],
            ["generate", "snakefile", "--help"],
            ["remote", "--help"],
            ["remote", "list", "--help"],
            ["remote", "add", "--help"],
            ["remote", "remove", "--help"],
            ["genome", "--help"],
            ["genome", "list", "--help"],
            ["genome", "remove", "--help"],
            ["serve", "--help"],
            ["dash", "--help"],
        ],
    )
    def test_help_exits_zero_with_output(self, cmd, capsys):
        from refgenie.cli.commands.framework import CleanHelpCliSource
        from refgenie.cli.parser import TopLevelParser

        with pytest.raises(SystemExit) as exc_info:
            cli_source = CleanHelpCliSource(TopLevelParser, cli_parse_args=cmd)
            TopLevelParser(_cli_settings_source=cli_source)
        assert exc_info.value.code == 0
        assert len(capsys.readouterr().out) >= 30, f"Help output too short for {cmd}"


# ===========================================================================
# Dispatch and exit codes (commands run against a Refgenie instance)
# ===========================================================================


# ---------------------------------------------------------------------------
# Friendly validation error messages
# ---------------------------------------------------------------------------


class TestValidationErrorMessages:
    """Missing/invalid required arguments produce friendly errors, not
    tracebacks or raw pydantic ValidationErrors."""

    def test_genome_init_missing_name(self, capsys):
        with pytest.raises(SystemExit) as exc_info:
            main(test_args=["genome", "init"])
        assert exc_info.value.code == 1
        err = capsys.readouterr().err
        assert "Missing required argument" in err
        assert "ValidationError" not in err
        assert "Traceback" not in err

    def test_help_hint_shown(self, capsys):
        with pytest.raises(SystemExit) as exc_info:
            main(test_args=["genome", "init"])
        assert exc_info.value.code == 1
        assert "--help" in capsys.readouterr().err

    def test_valid_command_no_validation_error(self, tmp_path, capsys):
        """A valid command must not raise a ValidationError (it may still exit
        for downstream reasons -- that's a different error)."""
        fasta = tmp_path / "test.fa"
        fasta.write_text(">chr1\nACGT\n")
        try:
            main(test_args=["genome", "init", "--name", "hg38", "--fasta", str(fasta)])
        except SystemExit:
            pass
        err = capsys.readouterr().err
        assert "Missing required argument" not in err
        assert "ValidationError" not in err

    def test_remote_add_no_args(self, capsys):
        with pytest.raises(SystemExit) as exc_info:
            main(test_args=["remote", "add"])
        assert exc_info.value.code == 1
        err = capsys.readouterr().err
        assert "ValidationError" not in err
        assert "Traceback" not in err
        assert "Missing required argument '--type'" in err
        assert "Missing required argument '--prefix'" in err
        assert "Missing required argument '--description'" in err
        assert "refgenie remote add --help" in err

    def test_remote_add_partial_args(self, capsys):
        with pytest.raises(SystemExit) as exc_info:
            main(test_args=["remote", "add", "--type", "s3"])
        assert exc_info.value.code == 1
        err = capsys.readouterr().err
        assert "Missing required argument '--prefix'" in err
        assert "Missing required argument '--description'" in err
        assert "Missing required argument '--type'" not in err

    def test_remote_add_invalid_type(self, capsys):
        with pytest.raises(SystemExit) as exc_info:
            main(test_args=["remote", "add", "--type", "ftp", "--prefix", "x", "--description", "y"])
        assert exc_info.value.code == 1
        err = capsys.readouterr().err
        assert "ValidationError" not in err
        assert "Traceback" not in err
        assert "--type" in err

    def test_remote_remove_missing_type(self, capsys):
        with pytest.raises(SystemExit) as exc_info:
            main(test_args=["remote", "remove"])
        assert exc_info.value.code == 1
        err = capsys.readouterr().err
        assert "Missing required argument '--type'" in err
        assert "refgenie remote remove --help" in err

    def test_rename_missing_new_name(self, capsys):
        with pytest.raises(SystemExit) as exc_info:
            main(test_args=["rename", "hg38/fasta"])
        assert exc_info.value.code == 1
        err = capsys.readouterr().err
        assert "Traceback" not in err
        # new_asset_name has alias 'n', so the flag should be '-n'
        assert "Missing required argument '-n'" in err


# ---------------------------------------------------------------------------
# Command dispatch on an empty (initialized) database
# ---------------------------------------------------------------------------


DC_ADD = ("data-channel", "add", "test_channel", "http", "http://example.com/index.yaml")
REMOTE_ADD = (
    "remote", "add", "--type", "http",
    "--prefix", "http://test.example.com", "--description", "test remote",
)


class TestEmptyDatabaseCommands:
    """Commands that work (or fail cleanly) without any built assets."""

    def test_init(self, cli, tmp_path):
        assert cli("init", "--genome-folder", str(tmp_path / "genomes")) == 0

    def test_list_empty(self, cli, cli_rg, capsys):
        assert cli("list") == 0
        err = capsys.readouterr().err
        assert err == "" or "No assets" in err

    def test_populate_with_file(self, cli, cli_rg, tmp_path):
        temp_file = tmp_path / "paths.txt"
        temp_file.write_text("rCRSd/fasta\n")
        assert cli("populate", "-f", str(temp_file)) == 0

    def test_populater_server_is_transient(self, cli, cli_rg, tmp_path):
        """`populater -s` is documented as not persisting; it must not subscribe."""
        temp_file = tmp_path / "paths.txt"
        temp_file.write_text("rCRSd/fasta\n")
        assert cli("populater", "-s", "http://transient.example.com", "-f", str(temp_file)) == 0
        assert list(cli_rg.sources.get_subscriptions()) == []

    @pytest.mark.parametrize(
        "argvs",
        [
            [("subscribe", "-s", "http://test-server.example.com")],
            [("subscribe", "-s", "http://a.com,http://b.com")],
            [("unsubscribe", "-s", "http://test-server.example.com")],
            [("purge", "--force")],
            [("recipe", "list")],
            [("asset-class", "list")],
            [("remote", "list")],
            [("genome", "list")],
            [REMOTE_ADD],
            [REMOTE_ADD, ("remote", "remove", "--type", "http")],
            [DC_ADD],
            [DC_ADD, ("data-channel", "remove", "test_channel")],
            [DC_ADD, ("data-channel", "show", "test_channel")],
        ],
        ids=[
            "subscribe", "subscribe-comma-list", "unsubscribe", "purge-force",
            "recipe-list", "asset-class-list", "remote-list", "genome-list",
            "remote-add", "remote-remove", "data-channel-add",
            "data-channel-remove", "data-channel-show",
        ],
    )
    def test_command_exits_zero(self, cli, cli_rg, argvs):
        """Each command exits 0; earlier entries in a row are setup."""
        for argv in argvs:
            assert cli(*argv) == 0, f"`refgenie {' '.join(argv)}` did not exit 0"

    @pytest.mark.parametrize(
        "argv",
        [
            ("config", "get"),
            ("alias", "get"),
            ("stage", "list"),
            ("data-channel", "list"),
        ],
        ids=["config-get", "alias-get", "stage-list", "data-channel-list"],
    )
    def test_command_exits_zero_with_output(self, cli, cli_rg, capsys, argv):
        assert cli(*argv) == 0
        assert len(capsys.readouterr().out) > 0

    @pytest.mark.parametrize(
        "argv",
        [
            ("recipe", "show", "nonexistent"),
            ("recipe", "remove", "nonexistent"),
            ("recipe", "requirements", "nonexistent"),
            ("asset-class", "show", "nonexistent"),
            ("asset-class", "remove", "nonexistent"),
            ("alias", "remove", "-a", "nonexistent_alias"),
            ("list", "-g", "g1,g2"),
            ("compare", "genome_a", "genome_b"),
        ],
        ids=[
            "recipe-show", "recipe-remove", "recipe-requirements",
            "asset-class-show", "asset-class-remove", "alias-remove",
            "list-unknown-genomes", "compare-unknown-genomes",
        ],
    )
    def test_nonexistent_target_exits(self, run_main, cli_rg, argv):
        with pytest.raises(SystemExit):
            run_main(*argv)

    @pytest.mark.parametrize(
        "argv",
        [
            ("config", "set"),
        ],
    )
    def test_unimplemented_command_exits(self, run_main, cli_rg, argv):
        with pytest.raises(SystemExit):
            run_main(*argv)


class TestArgumentValidationExitCodes:
    """Missing/invalid arguments produce a non-zero exit at dispatch time."""

    @pytest.mark.parametrize(
        "argv, code",
        [
            (("getseq", "-l", "chr1:1-10"), 1),
            (("getseq", "-g", "rCRSd"), 1),
            (("subscribe",), 1),
            (("data-channel", "add", "test", "invalid_type", "http://x.com"), 1),
            (("dash", "--port", "not_a_number"), 1),
            (("genome", "remove"), 1),
            (("compare", "genome1"), 2),
            (("add", "rCRSd/fasta:new", "-c", "fasta"), 1),
            (("add", "rCRSd/fasta:new", "-p", "/some/path"), 1),
        ],
        ids=[
            "getseq-no-genome", "getseq-no-locus", "subscribe-no-server",
            "data-channel-add-bad-type", "dash-bad-port", "genome-remove-no-arg",
            "compare-one-genome", "add-no-path", "add-no-asset-class",
        ],
    )
    def test_command_exits_nonzero(self, cli, cli_rg, argv, code):
        assert cli(*argv) == code


# ---------------------------------------------------------------------------
# Error contract: a failure must never exit 0
# ---------------------------------------------------------------------------


class TestPullExitCode:
    def test_failed_pull_exits_nonzero(self, cli, cli_rg):
        """`refgenie pull <nonexistent>` must not report success."""
        from refgenie.exceptions import PullFailedError

        with patch.object(
            Refgenie, "pull", side_effect=PullFailedError("asset not found on any server")
        ):
            assert cli("pull", "nonexistent/asset") != 0

    def test_pull_returning_none_exits_nonzero(self, cli, cli_rg):
        """Refgenie.pull signals failure by RETURNING None, not only by raising.

        This is the path a real `refgenie pull nonexistent/asset` takes: no
        server can serve the asset, AssetPuller logs and pull() returns None.
        """
        with patch.object(Refgenie, "pull", return_value=None):
            assert cli("pull", "nonexistent/asset") != 0

    def test_malformed_path_exits_nonzero(self, cli, cli_rg):
        assert cli("pull", "garbage") != 0

    @pytest.mark.parametrize(
        "argv",
        [
            ("recipe", "requirements", "fasta", "--recipe-version", "0.0.1"),
            ("build", "rCRSd/fasta", "--requirements", "--recipe-version", "0.0.1"),
        ],
        ids=["recipe-requirements", "build-requirements"],
    )
    def test_requirements_forwards_recipe_version(self, cli, cli_rg, argv):
        """Both --requirements paths must honor --recipe-version, not the latest."""
        from refgenie.managers.recipe import RecipeManager

        with patch.object(RecipeManager, "get", return_value=MagicMock(name="fasta")) as get:
            cli(*argv)
        assert get.call_args.kwargs["recipe_version"] == "0.0.1"

    def test_force_forwarded_on_registry_path(self, cli, cli_rg):
        """`pull <path> --force` must reach Refgenie.pull, not just --all/--asset."""
        with patch.object(Refgenie, "pull", return_value="ok") as mock_pull:
            cli("pull", "rCRSd/fasta", "--force")
        assert mock_pull.call_args.kwargs["force"] is True


class TestDataChannelHintOnError:
    """`pull` and `build` must explain an unsynced data channel, matching the
    hint `genome init` already prints (issue 12 of the migration testing
    report). The hint is byte-identical across all three call sites because
    they all go through `data_channel_hint()`."""

    def test_helper_message_contains_runnable_commands(self):
        from refgenie.cli.commands.helpers import data_channel_hint

        msg = data_channel_hint("refgenie pull rCRSd/fasta")
        assert "refgenie data_channel add refgenie https" in msg
        assert "refgenie data_channel sync refgenie --exists-ok" in msg
        assert msg.endswith("refgenie pull rCRSd/fasta")

    def test_helper_without_retry_command_omits_trailing_line(self):
        from refgenie.cli.commands.helpers import data_channel_hint

        msg = data_channel_hint()
        assert msg.count("\n") == 2

    def test_pull_missing_asset_class_on_empty_db_shows_hint(self, cli, cli_rg, caplog):
        """`pull` against a fresh, unsynced instance must show the same
        data-channel guidance genome init already gives."""
        from refgenie.exceptions import MissingAssetClassError

        with patch.object(Refgenie, "pull", side_effect=MissingAssetClassError("fasta")):
            with caplog.at_level("INFO"):
                assert cli("pull", "rCRSd/fasta") != 0
        assert "refgenie data_channel add" in caplog.text
        assert "refgenie data_channel sync" in caplog.text
        assert "refgenie pull rCRSd/fasta" in caplog.text

    def test_pull_missing_asset_class_with_channel_synced_omits_hint(
        self, cli, cli_rg, fixtures_path, caplog
    ):
        """When recipes/asset classes ARE registered, a genuinely missing one
        must not be misdirected to 'sync a data channel'."""
        from refgenie.exceptions import MissingAssetClassError

        register_fasta(cli_rg, fixtures_path)
        with patch.object(Refgenie, "pull", side_effect=MissingAssetClassError("bowtie2_index")):
            with caplog.at_level("INFO"):
                assert cli("pull", "rCRSd/bowtie2_index") != 0
        assert "refgenie data_channel add" not in caplog.text

    def test_build_missing_recipe_on_empty_db_shows_hint(self, cli, cli_rg, caplog):
        from refgenie.exceptions import MissingRecipeError

        with patch.object(Refgenie, "build_asset", side_effect=MissingRecipeError("fasta")):
            with caplog.at_level("INFO"):
                assert cli("build", "rCRSd/fasta") != 0
        assert "refgenie data_channel add" in caplog.text
        assert "refgenie data_channel sync" in caplog.text
        assert "refgenie build rCRSd/fasta" in caplog.text

    def test_build_missing_recipe_with_channel_synced_omits_hint(
        self, cli, cli_rg, fixtures_path, caplog
    ):
        from refgenie.exceptions import MissingRecipeError

        register_fasta(cli_rg, fixtures_path)
        with patch.object(Refgenie, "build_asset", side_effect=MissingRecipeError("custom_recipe")):
            with caplog.at_level("INFO"):
                assert cli("build", "rCRSd/custom_recipe") != 0
        assert "refgenie data_channel add" not in caplog.text


class TestSilentZeroPaths:
    """Handlers that used to log an error and fall off the end (exit 0)."""

    def test_data_channel_show_missing(self, cli, cli_rg):
        assert cli("data-channel", "show", "nosuchchannel") != 0

    def test_data_channel_add_failure(self, cli, cli_rg):
        # '__' in a channel name is rejected by SourceManager.add_channel.
        code = cli("data-channel", "add", "bad__name", "http", "http://example.com/index.yaml")
        assert code != 0

    def test_data_channel_remove_missing(self, cli, cli_rg):
        assert cli("data-channel", "remove", "nosuchchannel") != 0

    def test_data_channel_validate_invalid(self, cli, cli_rg):
        assert cli("data-channel", "validate", "nosuchchannel") != 0

    def test_remote_add_failure(self, cli, cli_rg):
        from refgenie.managers.configuration import ConfigurationManager

        with patch.object(ConfigurationManager, "add_remote", side_effect=ValueError("boom")):
            code = cli("remote", "add", "--type", "s3", "--prefix", "p", "--description", "d")
        assert code != 0

    def test_remote_remove_missing(self, cli, cli_rg):
        assert cli("remote", "remove", "--type", "s3") != 0

    def test_genome_sync_registration_failure(self, cli, cli_rg):
        """A sync where every collection fails to register must not exit 0."""
        import refgenie.cli.commands.genome as handlers_mod

        source = MagicMock()
        source.list_collections.return_value = {"results": [{"digest": "abc", "description": ""}]}
        with (
            patch.object(handlers_mod, "_make_source", return_value=source),
            patch.object(type(cli_rg.genome), "initialize_genome", side_effect=RuntimeError("nope")),
        ):
            code = cli("genome", "sync", "--server-url", "http://example.com")
        assert code != 0

    def test_genome_browse_all_servers_unreachable(self, cli, cli_rg):
        import refgenie.cli.commands.genome as handlers_mod

        with patch.object(handlers_mod, "_make_source", side_effect=ConnectionError("refused")):
            code = cli("genome", "browse", "--server-url", "http://example.com")
        assert code != 0


class TestCLIBulkPullIntegration:
    """The bulk pull CLI rejects incomplete argument combinations with exit code 1."""

    def test_pull_all_without_genome_fails(self, cli, cli_rg):
        assert cli("pull", "--all") == 1

    def test_pull_asset_without_genome_fails(self, cli, cli_rg):
        assert cli("pull", "--asset", "fasta") == 1

    def test_pull_all_genomes_all_blocked(self, cli, cli_rg):
        """--all-genomes --all is blocked (use mirror instead)."""
        assert cli("pull", "--all-genomes", "--all") == 1

    def test_pull_without_args_fails(self, cli, cli_rg):
        assert cli("pull") == 1


# ---------------------------------------------------------------------------
# Commands requiring a built asset (component tier)
# ---------------------------------------------------------------------------


@pytest.mark.component
class TestCommandsWithAssets:
    """CLI commands that require a built asset.

    Component tier: the autouse fixture builds a real FASTA asset on disk for
    every test (~0.23s apiece), combining real filesystem and SQLite.
    """

    @pytest.fixture(autouse=True)
    def setup_asset(self, cli_rg, fixtures_path):
        register_fasta(cli_rg, fixtures_path)
        cli_rg.genome.initialize_genome(
            fasta_file_path=fixtures_path / "rCRSd.fa",
            alias_names=["rCRSd"],
            description="rCRSd mitochondrial reference",
        )
        cli_rg.build_asset(
            recipe_name="fasta",
            genome_name="rCRSd",
            asset_group_name="fasta",
            asset_name="default",
        )

    def test_id_asset_digest_is_sha256(self, cli, capsys):
        """`id rCRSd/fasta` outputs a 64-char SHA-256 hex digest."""
        assert cli("id", "rCRSd/fasta") == 0
        digest = capsys.readouterr().out.strip()
        assert len(digest) == 64
        assert all(c in "0123456789abcdef" for c in digest)

    def test_id_genome_only(self, cli, capsys):
        """`id rCRSd` returns a 32-char genome digest."""
        assert cli("id", "rCRSd") == 0
        assert len(capsys.readouterr().out.strip()) == 32

    def test_build_requirements(self, cli, cli_rg_with_recipe, capsys):
        assert cli("build", "rCRSd/fasta", "--requirements") == 0
        assert "fasta" in capsys.readouterr().out.lower()

    def test_alias_set_no_genome_errors(self, cli):
        """alias set with a nonexistent digest exits non-zero, unless forced."""
        assert cli("alias", "set", "-a", "my_alias", "-d", "nonexistent_digest") != 0
        assert cli("alias", "set", "-a", "my_alias", "-d", "nonexistent_digest", "--force") == 0

    def test_remove_aliases_drops_alias_when_last_asset_goes(self, cli, cli_engine):
        """`remove --aliases` leaves no alias behind once nothing is left."""
        assert cli("remove", "rCRSd/fasta", "-f", "--aliases") == 0
        fresh = Refgenie(database_engine=cli_engine, suppress_migrations=True)
        assert [a.name for a in fresh.alias.list_all()] == []


# ---------------------------------------------------------------------------
# serve / dash / generate dispatch (mocked servers)
# ---------------------------------------------------------------------------


class TestServerDispatch:
    """serve and dash dispatch correctly without starting real servers."""

    def _dispatch(self, cli_engine, args, runner="run_server"):
        """Run main(args) with refgenie.server.main replaced by a mock;
        return the mock runner (`run_server` for serve, `run_local` for dash)."""
        import sys as _sys

        mock_runner = MagicMock()
        mock_module = MagicMock()
        setattr(mock_module, runner, mock_runner)
        with (
            patch.object(Refgenie, "get_default_database_engine", return_value=cli_engine),
            patch.dict(_sys.modules, {"refgenie.server.main": mock_module}),
            patch("webbrowser.open_new_tab", MagicMock()) as mock_open_tab,
        ):
            main(test_args=args)
        return mock_runner, mock_open_tab

    def test_serve_default(self, cli_engine, cli_rg):
        run_server, _ = self._dispatch(cli_engine, ["serve"])
        run_server.assert_called_once_with(port=8000, reload=False)

    def test_serve_custom_port(self, cli_engine, cli_rg):
        run_server, _ = self._dispatch(cli_engine, ["serve", "--port", "9999"])
        run_server.assert_called_once_with(port=9999, reload=False)

    def test_serve_reload_flag(self, cli_engine, cli_rg):
        run_server, _ = self._dispatch(cli_engine, ["serve", "--reload"])
        run_server.assert_called_once_with(port=8000, reload=True)

    def test_dash_default(self, cli_engine, cli_rg):
        """`dash` runs the local-mode app on the loopback address."""
        run_local, open_tab = self._dispatch(cli_engine, ["dash"], runner="run_local")
        run_local.assert_called_once_with(port=8080)
        open_tab.assert_called_once_with("http://127.0.0.1:8080")

    def test_dash_custom_port(self, cli_engine, cli_rg):
        run_local, open_tab = self._dispatch(
            cli_engine, ["dash", "--port", "7777"], runner="run_local"
        )
        run_local.assert_called_once_with(port=7777)
        open_tab.assert_called_once_with("http://127.0.0.1:7777")


class TestGenerateDispatch:
    def test_generate_snakefile(self, cli_engine, cli_rg, tmp_path):
        """generate snakefile dispatches to populate_snakefile_template with the
        requested output path."""
        mock_populate = MagicMock()
        output_path = tmp_path / "Snakefile"
        with (
            patch.object(Refgenie, "get_default_database_engine", return_value=cli_engine),
            patch("refgenie.snakefile.generate.populate_snakefile_template", mock_populate),
        ):
            main(test_args=["generate", "snakefile", "-o", str(output_path)])
        mock_populate.assert_called_once()
        call = mock_populate.call_args
        assert str(call.kwargs.get("snakefile_output_path")) == str(output_path)


# ---------------------------------------------------------------------------
# `refgenie id` and the remote-digest lookup behind it
#
# Covers the RefgetStore-enhanced ``handle_id`` features (--verbose,
# --validate-store, --info, --remote) and ``Refgenie.check_remote_digest``,
# which exists to serve ``id --remote``. Component tier: each test initializes
# a real genome from the rCRSd FASTA.
# ---------------------------------------------------------------------------


def _id_cmd(paths, *, verbose=False, validate_store=False, remote=False, info=False):
    """Build a MagicMock command namespace for handle_id."""
    cmd = MagicMock()
    cmd.asset_registry_paths = paths
    cmd.verbose = verbose
    cmd.validate_store = validate_store
    cmd.remote = remote
    cmd.info = info
    return cmd


def _capture_handle_id(cmd, rgc):
    """Run handle_id and capture stdout (handles rich.print redirects)."""
    from refgenie.cli.commands.lookup import handle_id

    captured = io.StringIO()
    old_stdout = sys.stdout
    sys.stdout = captured
    try:
        handle_id(cmd, rgc)
    finally:
        sys.stdout = old_stdout
    return captured.getvalue()


@pytest.mark.component
class TestIdVerboseAndValidate:
    def test_verbose_output_contains_metadata(self, refgenie_with_genome):
        cmd = _id_cmd(["rCRSd"], verbose=True)
        output = _capture_handle_id(cmd, refgenie_with_genome)
        assert "digest:" in output
        assert "sequences:" in output
        assert "total_length:" in output
        assert "source:" in output

    def test_validate_succeeds_for_valid_genome(self, refgenie_with_genome):
        cmd = _id_cmd(["rCRSd"], validate_store=True)
        output = _capture_handle_id(cmd, refgenie_with_genome)
        assert len(output.strip()) > 0


@pytest.mark.component
class TestIdInfo:
    def test_info_shows_aliases_and_metadata(self, refgenie_with_genome):
        digest = refgenie_with_genome.alias.resolve("rCRSd")
        output = _capture_handle_id(_id_cmd([digest], info=True), refgenie_with_genome)
        assert "digest:" in output
        assert "aliases:" in output
        assert "rCRSd" in output
        assert "sequences:" in output
        assert "total_length:" in output

    def test_info_unknown_digest_errors(self, refgenie_with_genome):
        cmd = _id_cmd(["nonexistent_digest_0000000000000"], info=True)
        with pytest.raises(SystemExit):
            _capture_handle_id(cmd, refgenie_with_genome)


@pytest.mark.component
class TestIdRemote:
    def test_remote_with_unresolvable_name_errors(self, refgenie_with_genome):
        cmd = _id_cmd(["GRCh38"], remote=True)  # not a known alias
        with pytest.raises(SystemExit):
            _capture_handle_id(cmd, refgenie_with_genome)

    def test_remote_with_digest_calls_check_remote(self, refgenie_with_genome):
        rgc = refgenie_with_genome
        fake_digest = "a" * 32
        cmd = _id_cmd([fake_digest], remote=True)
        with patch.object(
            rgc, "check_remote_digest", return_value={"digest": fake_digest}
        ) as mock_check:
            output = _capture_handle_id(cmd, rgc)
        mock_check.assert_called_once_with(fake_digest)
        assert fake_digest in output

    def test_remote_with_digest_not_found_errors(self, refgenie_with_genome):
        rgc = refgenie_with_genome
        cmd = _id_cmd(["b" * 32], remote=True)
        with patch.object(rgc, "check_remote_digest", return_value=None):
            with pytest.raises(SystemExit):
                _capture_handle_id(cmd, rgc)


@pytest.mark.component
class TestCheckRemoteDigest:
    def test_returns_none_with_no_servers(self, refgenie_with_genome):
        assert refgenie_with_genome.check_remote_digest("some_digest", remote_servers=[]) is None

    def test_returns_none_on_exception(self, refgenie_with_genome):
        with patch(
            "refgenie.core.sequences.make_source",
            side_effect=ConnectionError("unreachable"),
        ) as mock_make_source:
            result = refgenie_with_genome.check_remote_digest(
                "some_digest", remote_servers=["http://fake.server"]
            )
        assert result is None
        mock_make_source.assert_called_once_with("http://fake.server")

    def test_returns_none_when_source_lacks_collection(self, refgenie_with_genome):
        source = MagicMock()
        source.verify_collection.return_value = None
        with patch("refgenie.core.sequences.make_source", return_value=source):
            result = refgenie_with_genome.check_remote_digest(
                "some_digest", remote_servers=["http://fake.server"]
            )
        assert result is None
        source.verify_collection.assert_called_once_with("some_digest")


# ---------------------------------------------------------------------------
# `refgenie-build-fasta` console script (refgenie.cli.build_fasta)
#
# This is a separate entry point from the main ``refgenie`` CLI: it renders a
# FASTA, its .fai index and a .chrom.sizes file out of a RefgetStore. Unit
# tier -- the store is built in tmp_path and no external binary is involved.
# ---------------------------------------------------------------------------


@pytest.fixture
def store_with_genome(tmp_path, fixtures_path):
    """A RefgetStore populated with the rCRSd test FASTA."""
    from gtars.refget import RefgetStore

    store_path = tmp_path / "store"
    store = RefgetStore.on_disk(store_path)
    metadata, _was_new = store.add_sequence_collection_from_fasta(fixtures_path / "rCRSd.fa")
    return store_path, metadata.digest


class TestBuildFastaCli:
    def test_produces_all_files(self, store_with_genome, tmp_path):
        from refgenie.cli.build_fasta import build_fasta_main

        store_path, digest = store_with_genome
        output = tmp_path / "output"
        build_fasta_main([str(store_path), digest, str(output)])
        assert (output / f"{digest}.fa").exists()
        assert (output / f"{digest}.fa.fai").exists()
        assert (output / f"{digest}.chrom.sizes").exists()

    def test_fa_contains_sequences(self, store_with_genome, tmp_path):
        from refgenie.cli.build_fasta import build_fasta_main

        store_path, digest = store_with_genome
        output = tmp_path / "output"
        build_fasta_main([str(store_path), digest, str(output)])
        fa_content = (output / f"{digest}.fa").read_text()
        assert fa_content.startswith(">")
        assert len(fa_content) > 100

    def test_fai_has_five_columns(self, store_with_genome, tmp_path):
        from refgenie.cli.build_fasta import build_fasta_main

        store_path, digest = store_with_genome
        output = tmp_path / "output"
        build_fasta_main([str(store_path), digest, str(output)])
        fai_lines = (output / f"{digest}.fa.fai").read_text().strip().split("\n")
        assert fai_lines
        for line in fai_lines:
            assert len(line.split("\t")) == 5

    def test_chrom_sizes_matches_level2(self, store_with_genome, tmp_path):
        from gtars.refget import RefgetStore
        from refgenie.cli.build_fasta import build_fasta_main

        store_path, digest = store_with_genome
        output = tmp_path / "output"
        build_fasta_main([str(store_path), digest, str(output)])
        level2 = RefgetStore.open_local(store_path).get_collection_level2(digest)
        lines = (output / f"{digest}.chrom.sizes").read_text().strip().split("\n")
        assert lines
        for i, line in enumerate(lines):
            name, length = line.split("\t")
            assert name == level2["names"][i]
            assert int(length) == level2["lengths"][i]

    def test_custom_line_width(self, store_with_genome, tmp_path):
        """--line-width wraps FASTA sequence lines at the requested width."""
        from refgenie.cli.build_fasta import build_fasta_main

        store_path, digest = store_with_genome
        output = tmp_path / "output"
        build_fasta_main([str(store_path), digest, str(output), "--line-width", "60"])
        seq_lines = [
            line
            for line in (output / f"{digest}.fa").read_text().strip().split("\n")
            if not line.startswith(">")
        ]
        # rCRSd is ~33 kb, so wrapping at 60 yields many full-width lines.
        assert len(seq_lines) > 1
        assert len(seq_lines[0]) == 60

    def test_invalid_digest_fails(self, store_with_genome, tmp_path):
        from refgenie.cli.build_fasta import build_fasta_main

        store_path, _ = store_with_genome
        output = tmp_path / "output"
        # gtars surfaces an unknown collection digest as OSError from
        # RefgetStore.load_collection.
        with pytest.raises(OSError, match="Error loading collection"):
            build_fasta_main([str(store_path), "nonexistent_digest", str(output)])
