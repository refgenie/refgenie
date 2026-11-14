"""The refgenie CLI: argv parsing, --help output, dispatch and exit codes."""

import inspect
import io
import re
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from pydantic import ValidationError
from pydantic_settings import get_subcommand

from refgenie import Refgenie
from refgenie.cli import main as cli_main
from refgenie.cli.main import main
from refgenie.cli.commands import asset as cli_asset
from refgenie.cli.commands import helpers as cli_helpers
from refgenie.cli.commands import listing as cli_listing
from refgenie.cli.commands.pull import MirrorModel, PullModel
from refgenie.cli.dispatch import get_dispatch
from refgenie.managers.build import BuildManager
from refgenie.cli.errors import EXIT_INVALID_INPUT, EXIT_NOT_FOUND
from refgenie.managers.sources.servers import ServerManager
from refgenie.managers.transfer import TransferManager
from refgenie.managers.transfer.puller import BulkPullResult
from refgenie.models import GenomeAlias, GenomeDigest
from tests.support.cli import parse, subcmd
from tests.support.factories import fake_digest, make_engine, register_fasta


# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------


@pytest.fixture
def cli_engine():
    """In-memory SQLite engine for CLI tests."""
    return make_engine()


@pytest.fixture
def cli_rg(cli_engine, tmp_path):
    """Initialized Refgenie instance sharing the CLI engine."""
    rg = Refgenie(database_engine=cli_engine, suppress_migrations=True)
    rg.database.init(genome_folder=tmp_path / "genomes")
    return rg


@pytest.fixture
def cli_rg_with_recipe(tmp_path):
    """Initialized Refgenie with its own engine (used by build --requirements)."""
    rg = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    rg.database.init(genome_folder=tmp_path / "genomes")
    return rg


@pytest.fixture
def cli(cli_engine):
    """CLI invoker: invoke(*args) -> exit_code (use capsys to capture output)."""

    def invoke(*args):
        with patch("refgenie.core.root.create_database_engine", return_value=cli_engine):
            try:
                main(test_args=list(args))
                return 0
            except SystemExit as e:
                return e.code or 0

    return invoke


# ===========================================================================
# argv -> pydantic parsing and model validation (no database, no filesystem)
# ===========================================================================


class TestSubcommandParsing:
    def test_top_level_subcommands_parse(self):
        """Each top-level subcommand routes to its model with the right fields."""
        assert parse(["seek", "hg38/fasta"]).seek.asset_registry_paths == ["hg38/fasta"]

        cmp = parse(["compare", "hg38", "mm10"]).compare
        assert cmp.genomes == ["hg38", "mm10"]
        cmp = parse(["compare", "hg38", "--genome-digest", "d" * 32]).compare
        assert cmp.genomes == ["hg38"] and cmp.genome_digest == ["d" * 32]

        rm = parse(["remove", "hg38/fasta", "-f"]).remove
        assert rm.force is True and rm.asset_registry_paths == ["hg38/fasta"]

        assert parse(["listr", "-s", "http://example.com"]).listr.genome_server == [
            "http://example.com"
        ]
        assert parse(["id", "hg38"]).id.asset_registry_paths == ["hg38"]
        assert parse(["rename", "hg38/fasta", "-n", "fasta_v2"]).rename.new_asset_name == "fasta_v2"

        gs = parse(["getseq", "-g", "hg38", "-l", "chr1:1-100"]).getseq
        assert gs.genome == "hg38" and gs.locus == "chr1:1-100"

        add = parse(["add", "hg38/fasta:custom", "-p", "/tmp/path", "-c", "fasta"]).add
        assert add.asset_registry_paths == ["hg38/fasta:custom"]
        assert add.path == "/tmp/path" and add.asset_class == "fasta"

        assert parse(["populate", "-f", "/tmp/file.txt"]).populate.file == "/tmp/file.txt"

        seekr = parse(["seekr", "hg38/fasta", "-s", "http://example.com"]).seekr
        assert seekr.asset_registry_paths == ["hg38/fasta"]
        assert seekr.genome_server == ["http://example.com"]

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
        assert parse(["asset-class", "show", "fasta"]).asset_class.show.asset_class_name == "fasta"
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
            ["remote", "add", "--type", "http", "--prefix", "http://x.com", "--name", "test"]
        ).remote.add
        assert r_add.type == "http" and r_add.prefix == "http://x.com"
        assert r_add.name == "test"

        assert parse(["remote", "remove", "test"]).remote.remove.remote == "test"

        assert parse(["genome", "list"]).genome.list is not None

        g_rm = parse(["genome", "remove", "--genome", "rCRSd", "-f"]).genome.remove
        assert g_rm.genome == ["rCRSd"] and g_rm.force is True

    def test_build_container_flag(self):
        assert subcmd(["build", "rCRSd/fasta", "--container"]).container is True

    @pytest.mark.parametrize(
        "argv, expected",
        [
            (["--files", "fasta=/a.fa"], ["fasta=/a.fa"]),
            (["--files", "fasta=/a.fa", "--files", "gtf=/b.gtf"], ["fasta=/a.fa", "gtf=/b.gtf"]),
        ],
        ids=["one", "repeated"],
    )
    def test_build_files_flag(self, argv, expected):
        assert subcmd(["build", "rCRSd/fasta", *argv]).files == expected


class TestMutuallyExclusiveFlags:
    @pytest.mark.parametrize(
        "argv, match",
        [
            (["pull", "--skip-large", "--pull-large"], "mutually exclusive"),
            (["pull", "--batch", "--skip-large"], "conflicts with --skip-large"),
            (["mirror", "--skip-large", "--pull-large"], "mutually exclusive"),
            (["data-channel", "sync", "ch", "--exists-ok", "--exists-overwrite"], "mutually"),
            (["alias", "get", "-a", "x", "-g", "y"], "mutually exclusive"),
        ],
        ids=["pull-large", "pull-batch", "mirror-large", "data-channel-sync", "alias-get"],
    )
    def test_conflicting_flags_are_rejected(self, argv, match):
        with pytest.raises(ValidationError, match=match):
            parse(argv)

    @pytest.mark.parametrize(
        "argv, force_large",
        [([], None), (["--skip-large"], False), (["--pull-large"], True), (["--batch"], True)],
        ids=["default", "skip-large", "pull-large", "batch"],
    )
    def test_large_archive_choice(self, argv, force_large):
        assert subcmd(["pull", *argv]).resolve_force_large() is force_large


class TestPullArgParsing:
    @pytest.mark.parametrize(
        "argv, model, expected",
        [
            (["pull", "GRCm39/fasta"], PullModel, {"asset_registry_paths": ["GRCm39/fasta"]}),
            (
                ["pull", "GRCm39/fasta", "GRCh38/fasta"],
                PullModel,
                {"asset_registry_paths": ["GRCm39/fasta", "GRCh38/fasta"]},
            ),
            (
                ["pull", "GRCm39/fasta", "--force"],
                PullModel,
                {"asset_registry_paths": ["GRCm39/fasta"], "force": True},
            ),
            (
                ["pull", "-g", "hg38,mm10", "--all"],
                PullModel,
                {"genome": ["hg38", "mm10"], "all": True, "asset_registry_paths": []},
            ),
            (
                ["pull", "--all-genomes", "--asset", "fasta"],
                PullModel,
                {"all_genomes": True, "asset": "fasta"},
            ),
            (["pull", "-g", "hg38", "--init"], PullModel, {"genome": ["hg38"], "init": True}),
            (
                ["pull", "--genome-digest", "d" * 32, "--init"],
                PullModel,
                {"genome_digest": ["d" * 32], "init": True},
            ),
            (["pull", "-g", "hg38", "--all", "--force"], PullModel, {"force": True}),
            (["pull", "--size-cutoff", "5.5"], PullModel, {"size_cutoff": 5.5}),
            (
                ["pull", "fasta", "--genome-digest", "D"],
                PullModel,
                {"asset_registry_paths": ["fasta"], "genome_digest": ["D"]},
            ),
            (["pull", "--batch"], PullModel, {"batch": True, "pull_large": True}),
            (["mirror"], MirrorModel, {}),
            (["mirror", "--force"], MirrorModel, {"force": True}),
            (["mirror", "--size-cutoff", "20"], MirrorModel, {"size_cutoff": 20.0}),
        ],
        ids=[
            "pull-positional-path",
            "pull-positional-paths",
            "pull-positional-force",
            "pull-all-comma-genomes",
            "pull-all-genomes-asset",
            "pull-init",
            "pull-init-digest",
            "pull-all-force",
            "pull-size-cutoff",
            "pull-positional-digest",
            "pull-batch",
            "mirror",
            "mirror-force",
            "mirror-size-cutoff",
        ],
    )
    def test_pull_arg_parsing(self, argv, model, expected):
        parsed = subcmd(argv)
        assert isinstance(parsed, model)
        for name, value in expected.items():
            assert getattr(parsed, name) == value, f"{name} mismatch for {' '.join(argv)}"


class TestAssetListAliasesTopLevelList:
    """`refgenie asset list` forwards to handle_list, so it must deserialize
    `-g a,b` the same way `refgenie list` does."""

    def test_comma_separated_genomes_match(self):
        top = parse(["list", "-g", "a,b"]).list
        nested = parse(["asset", "list", "-g", "a,b"]).asset.list
        assert nested.genome == top.genome == ["a", "b"]

    def test_genome_field_deserializes_like_list_model(self):
        """The two models must accept the same inputs, not just the same argv.

        Both fields must be CliList: a plain list[str] works through argparse
        but rejects a comma-separated string on direct construction.
        """
        from refgenie.cli.commands.asset import AssetListNestedModel
        from refgenie.cli.commands.listing import ListModel

        assert AssetListNestedModel(genome="a,b").genome == ListModel(genome="a,b").genome

    def test_asset_list_subcommand_routes_to_nested_model(self):
        """``refgenie asset list`` parses to AssetListNestedModel."""
        from refgenie.cli.commands.asset import AssetListNestedModel

        asset_group = subcmd(["asset", "list"])
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


@pytest.fixture(scope="module")
def help_output():
    buf = io.StringIO()
    with redirect_stdout(buf):
        with pytest.raises(SystemExit):
            main(["--help"])
    return buf.getvalue()


class TestTopLevelHelp:
    """`refgenie --help` groups commands under headers from COMMAND_GROUPS."""

    def test_all_commands_present(self, help_output):
        from refgenie.cli.messages import COMMAND_GROUPS

        for commands in COMMAND_GROUPS.values():
            for cmd in commands:
                assert cmd in help_output, f"Command '{cmd}' missing from --help"

    def test_group_headers_before_commands(self, help_output):
        """Each group header appears before all its member commands in the
        grouped section (searched after the usage/options block)."""
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


class TestPullHelp:
    def test_pull_help_shows_the_positional(self, capsys):
        with pytest.raises(SystemExit):
            main(["pull", "--help"])
        out = capsys.readouterr().out
        assert "ASSET-REGISTRY-PATHS" in out
        assert "--asset-registry-paths" not in out


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
        from refgenie.cli.framework import CleanHelpCliSource
        from refgenie.cli.parser import TopLevelParser

        with pytest.raises(SystemExit) as exc_info:
            cli_source = CleanHelpCliSource(TopLevelParser, cli_parse_args=cmd)
            TopLevelParser(_cli_settings_source=cli_source)
        assert exc_info.value.code == 0
        assert len(capsys.readouterr().out) >= 30, f"Help output too short for {cmd}"


# ===========================================================================
# Dispatch and exit codes (commands run against a Refgenie instance)
# ===========================================================================


class TestErrorExitCodes:
    """Bad input exits with the right code and a readable message, never a traceback."""

    @pytest.mark.parametrize(
        "argv, code, fragments",
        [
            # Parse errors: missing or invalid arguments.
            (("genome", "init"), 1, ["Missing required argument", "--help"]),
            (
                ("remote", "add"),
                1,
                [
                    "Missing required argument '--type'",
                    "Missing required argument '--prefix'",
                    "Missing required argument '--name'",
                    "refgenie remote add --help",
                ],
            ),
            (("remote", "add", "--type", "ftp", "--prefix", "x", "--name", "y"), 1, ["--type"]),
            (("remote", "remove"), 2, ["the following arguments are required: REMOTE"]),
            (("rename", "hg38/fasta"), 1, ["Missing required argument '-n'"]),
            (("dash", "--port", "not_a_number"), 1, []),
            (("data-channel", "add", "test", "invalid_type", "http://x.com"), 1, []),
            # Missing targets.
            (("recipe", "show", "nonexistent"), 1, ["Recipe 'nonexistent' not found"]),
            (("recipe", "remove", "nonexistent"), 1, ["Recipe 'nonexistent' not found"]),
            (("recipe", "requirements", "nonexistent"), 1, ["Recipe 'nonexistent' not found"]),
            (("asset-class", "show", "nonexistent"), 1, ["Asset class 'nonexistent' not found"]),
            (("asset-class", "remove", "nonexistent"), 1, ["Asset class 'nonexistent'"]),
            (("alias", "remove", "-a", "nonexistent_alias"), 1, ["'nonexistent_alias' not found"]),
            (("list", "-g", "g1,g2"), 1, ["Genome alias 'g1' not found"]),
            (
                ("compare", "genome_a", "genome_b"),
                EXIT_NOT_FOUND,
                ["Genome 'genome_a' is not registered here"],
            ),
            (("data-channel", "show", "nosuchchannel"), 2, ["'nosuchchannel' not found"]),
            (("data-channel", "remove", "nosuchchannel"), 2, ["'nosuchchannel' not found"]),
            (("data-channel", "validate", "nosuchchannel"), 1, ["'nosuchchannel' is invalid"]),
            (("data-channel", "add", "bad__name", "http", "http://x.com/i.yaml"), 1, ["'__'"]),
            (("remote", "remove", "no-such-remote"), 2, ["'no-such-remote' not found"]),
            # Incomplete argument combinations.
            (("getseq", "-l", "chr1:1-10"), 3, []),
            (("getseq", "-g", "rCRSd"), 1, []),
            (("subscribe",), 1, []),
            (("genome", "remove"), 3, []),
            (("compare", "genome1"), 3, []),
            (("add", "rCRSd/fasta:new", "-c", "fasta"), 1, []),
            (("add", "rCRSd/fasta:new", "-p", "/some/path"), 1, []),
            (
                ("add", "fasta:new", "-c", "fasta", "-p", "/some/path"),
                EXIT_INVALID_INPUT,
                ["genome/asset_group[.seek_key][:asset]"],
            ),
            (("pull",), 1, ["Specify assets to pull"]),
            (("pull", "--all"), 1, ["Specify genome(s)"]),
            (("pull", "--asset", "fasta"), 1, ["Specify genome(s)"]),
            (("pull", "--all-genomes", "--all"), 1, ["Use 'refgenie mirror'"]),
            (("pull", "garbage"), EXIT_INVALID_INPUT, ["Invalid asset path: 'garbage'"]),
        ],
        ids=lambda v: "-".join(v) if isinstance(v, tuple) else None,
    )
    def test_exit_code_and_message(self, cli, cli_rg, capsys, caplog, argv, code, fragments):
        assert cli(*argv) == code
        text = capsys.readouterr().err + caplog.text
        assert "Traceback" not in text and "ValidationError" not in text
        for fragment in fragments:
            assert fragment in text

    def test_only_the_missing_arguments_are_named(self, cli, capsys):
        assert cli("remote", "add", "--type", "s3") == 1
        err = capsys.readouterr().err
        assert "Missing required argument '--prefix'" in err
        assert "Missing required argument '--type'" not in err


DC_ADD = ("data-channel", "add", "test_channel", "http", "http://example.com/index.yaml")
REMOTE_ADD = (
    "remote",
    "add",
    "--type",
    "http",
    "--prefix",
    "http://test.example.com",
    "--name",
    "test-remote",
)


class TestEmptyDatabaseCommands:
    """Commands that work without any built assets."""

    def test_init(self, cli, tmp_path):
        assert cli("init", "--genome-folder", str(tmp_path / "genomes")) == 0

    def test_init_subscribes_default_server(self, cli, cli_engine, tmp_path, monkeypatch):
        from refgenie.config import config
        from refgenie.const import DEFAULT_SERVER_URL

        monkeypatch.setattr(config, "subscribe_default_server", True)
        assert cli("init", "--genome-folder", str(tmp_path / "genomes")) == 0
        rg = Refgenie(database_engine=cli_engine, suppress_migrations=True)
        assert rg.servers.subscriptions() == [DEFAULT_SERVER_URL]

    def test_init_no_default_server(self, cli, cli_engine, tmp_path, monkeypatch):
        from refgenie.config import config

        monkeypatch.setattr(config, "subscribe_default_server", True)
        assert cli("init", "--genome-folder", str(tmp_path / "g"), "--no-default-server") == 0
        rg = Refgenie(database_engine=cli_engine, suppress_migrations=True)
        assert rg.servers.subscriptions() == []

    def test_list_empty(self, cli, cli_rg, capsys):
        assert cli("list") == 0
        err = capsys.readouterr().err
        assert err == "" or "No assets" in err

    def test_seek_missing_asset_logs_to_stderr_only(self, cli, cli_rg, capsys):
        """Scripts read `seek` stdout as a path; errors must not land there."""
        assert cli("seek", "nosuchgenome/fasta") == EXIT_NOT_FOUND
        captured = capsys.readouterr()
        assert captured.out == ""
        assert "nosuchgenome" in captured.err

    def test_populate_with_file(self, cli, cli_rg, tmp_path):
        temp_file = tmp_path / "paths.txt"
        temp_file.write_text("rCRSd/fasta\n")
        assert cli("populate", "-f", str(temp_file)) == 0

    def test_populater_server_is_transient(self, cli, cli_rg, tmp_path):
        """`populater -s` is documented as not persisting; it must not subscribe."""
        temp_file = tmp_path / "paths.txt"
        temp_file.write_text("rCRSd/fasta\n")
        assert cli("populater", "-s", "http://transient.example.com", "-f", str(temp_file)) == 0
        assert cli_rg.servers.subscriptions() == []

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
            [REMOTE_ADD, ("remote", "remove", "test-remote")],
            [REMOTE_ADD, ("remote", "remove", "1")],
            [DC_ADD],
            [DC_ADD, ("data-channel", "remove", "test_channel")],
            [DC_ADD, ("data-channel", "show", "test_channel")],
        ],
        ids=[
            "subscribe",
            "subscribe-comma-list",
            "unsubscribe",
            "purge-force",
            "recipe-list",
            "asset-class-list",
            "remote-list",
            "genome-list",
            "remote-add",
            "remote-remove-by-name",
            "remote-remove-by-id",
            "data-channel-add",
            "data-channel-remove",
            "data-channel-show",
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


# ---------------------------------------------------------------------------
# Error contract: a failure must never exit 0
# ---------------------------------------------------------------------------


class TestPullExitCode:
    def test_failed_pull_exits_nonzero(self, cli, cli_rg):
        """`refgenie pull <nonexistent>` must not report success."""
        from refgenie.exceptions import PullFailedError

        with patch.object(
            TransferManager, "pull", side_effect=PullFailedError("asset not found on any server")
        ):
            assert cli("pull", "nonexistent/asset") != 0

    def test_pull_returning_none_exits_nonzero(self, cli, cli_rg):
        """transfer.pull returns None for a skipped pull; nothing was retrieved,
        so that is not a success either."""
        with patch.object(TransferManager, "pull", return_value=None):
            assert cli("pull", "nonexistent/asset") != 0

    def test_server_error_after_alias_resolution_exits_nonzero(self, cli, cli_rg, caplog):
        """The genome resolves, the subscribed server then fails (a 500 from its
        relationships endpoint): the pull exits non-zero and says what the
        server said, not "no subscribed servers"."""
        from refgenie.managers.sources.servers import ServerManager
        from refgenie.managers.transfer.puller import AssetPuller

        url = "http://broken-server"
        cli_rg.servers.subscribe([url])
        digest = fake_digest("g")
        client = MagicMock()
        client.asset_group.list.return_value = [
            {"id": "group1", "name": "bowtie2_index", "genome_digest": digest}
        ]
        client.asset.list.return_value = [
            {"id": "asset1", "name": "default", "digest": "d1", "serving_modes": ["archive"]}
        ]
        client.asset.relationships.side_effect = RuntimeError(
            "500 Server Error: Internal Server Error"
        )

        with (
            patch.object(
                AssetPuller,
                "_resolve_pull_genome",
                return_value=(GenomeAlias("hg38"), GenomeDigest(digest), {}),
            ),
            patch.object(ServerManager, "client", return_value=client),
            caplog.at_level("INFO"),
        ):
            assert cli("pull", "hg38/bowtie2_index") == EXIT_NOT_FOUND

        assert "no server could serve 'hg38/bowtie2_index'" in caplog.text
        assert f"{url}: Failed to query asset relationships" in caplog.text
        assert "500 Server Error" in caplog.text
        assert "No server subscriptions" not in caplog.text

    def test_no_subscriptions_exits_nonzero_with_hint(self, cli, cli_rg, caplog):
        """Nothing subscribed and nobody to say yes to the default server."""
        with caplog.at_level("INFO"):
            assert cli("pull", "hg38/fasta") == EXIT_NOT_FOUND
        assert "refgenie subscribe" in caplog.text

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
        """`pull <path> --force` must reach transfer.pull, not just --all/--asset."""
        with patch.object(TransferManager, "pull", return_value="ok") as mock_pull:
            cli("pull", "rCRSd/fasta", "--force")
        assert mock_pull.call_args.kwargs["force"] is True


class TestDataChannelHintOnError:
    """`pull` and `build` explain an unsynced data channel with the same hint
    `genome init` prints, via `data_channel_hint()`."""

    def test_helper_message_contains_runnable_commands(self):
        from refgenie.cli.commands.helpers import data_channel_hint

        msg = data_channel_hint("refgenie pull rCRSd/fasta")
        assert "refgenie data-channel add refgenie https" in msg
        assert "refgenie data-channel sync refgenie --exists-ok" in msg
        assert msg.endswith("refgenie pull rCRSd/fasta")

    def test_helper_without_retry_command_omits_trailing_line(self):
        from refgenie.cli.commands.helpers import data_channel_hint

        msg = data_channel_hint()
        assert msg.count("\n") == 2

    def test_pull_missing_asset_class_on_empty_db_shows_hint(self, cli, cli_rg, capsys):
        """`pull` against a fresh, unsynced instance must show the same
        data-channel guidance genome init already gives."""
        from refgenie.exceptions import MissingAssetClassError

        with patch.object(TransferManager, "pull", side_effect=MissingAssetClassError("fasta")):
            assert cli("pull", "rCRSd/fasta") != 0
        err = capsys.readouterr().err
        assert "refgenie data-channel add" in err
        assert "refgenie data-channel sync" in err
        assert "refgenie pull rCRSd/fasta" in err

    def test_pull_missing_asset_class_with_channel_synced_omits_hint(
        self, cli, cli_rg, fixtures_path, capsys
    ):
        """When recipes/asset classes ARE registered, a genuinely missing one
        must not be misdirected to 'sync a data channel'."""
        from refgenie.exceptions import MissingAssetClassError

        register_fasta(cli_rg, fixtures_path)
        with patch.object(
            TransferManager, "pull", side_effect=MissingAssetClassError("bowtie2_index")
        ):
            assert cli("pull", "rCRSd/bowtie2_index") != 0
        err = capsys.readouterr().err
        assert "refgenie data-channel add" not in err

    def test_build_missing_recipe_on_empty_db_shows_hint(self, cli, cli_rg, capsys):
        from refgenie.exceptions import MissingRecipeError

        with patch.object(BuildManager, "run", side_effect=MissingRecipeError("fasta")):
            assert cli("build", "rCRSd/fasta") != 0
        err = capsys.readouterr().err
        assert "refgenie data-channel add" in err
        assert "refgenie data-channel sync" in err
        assert "refgenie build rCRSd/fasta" in err

    def test_build_missing_recipe_with_channel_synced_omits_hint(
        self, cli, cli_rg, fixtures_path, capsys
    ):
        from refgenie.exceptions import MissingRecipeError

        register_fasta(cli_rg, fixtures_path)
        with patch.object(BuildManager, "run", side_effect=MissingRecipeError("custom_recipe")):
            assert cli("build", "rCRSd/custom_recipe") != 0
        err = capsys.readouterr().err
        assert "refgenie data-channel add" not in err


class TestSilentZeroPaths:
    """A handler that logs an error must also exit nonzero, not fall off the end
    (exit 0)."""

    def test_build_failure(self, cli, cli_rg, fixtures_path):
        register_fasta(cli_rg, fixtures_path)
        with patch.object(BuildManager, "run", return_value=None):
            assert cli("build", "rCRSd/fasta") != 0

    def test_remote_add_failure(self, cli, cli_rg):
        from refgenie.managers.remote import RemoteManager

        with patch.object(RemoteManager, "add", side_effect=ValueError("boom")):
            code = cli("remote", "add", "--type", "s3", "--prefix", "p", "--name", "d")
        assert code != 0

    def test_genome_sync_registration_failure(self, cli, cli_rg):
        """A sync where every collection fails to register must not exit 0."""
        import refgenie.cli.commands.genome as handlers_mod

        source = MagicMock()
        source.list_collections.return_value = {"results": [{"digest": "abc", "description": ""}]}
        with (
            patch.object(handlers_mod, "make_source", return_value=source),
            patch.object(
                type(cli_rg.genome), "initialize_genome", side_effect=RuntimeError("nope")
            ),
        ):
            code = cli("genome", "sync", "--server-url", "http://example.com")
        assert code != 0

    def test_genome_browse_all_servers_unreachable(self, cli, cli_rg):
        import refgenie.cli.commands.genome as handlers_mod

        with patch.object(handlers_mod, "make_source", side_effect=ConnectionError("refused")):
            code = cli("genome", "browse", "--server-url", "http://example.com")
        assert code != 0


class TestBulkPullExitCode:
    @pytest.mark.parametrize("mode", [["--all"], ["--asset", "fasta"]])
    def test_bulk_pull_that_resolves_nothing_fails(self, cli, cli_rg, mode, caplog):
        """A bulk pull whose genome no server knows exits nonzero, not 0."""
        with patch.object(type(cli_rg.genome), "init_from_remote", return_value=None):
            assert cli("pull", "-g", "nosuchgenome", *mode) == EXIT_NOT_FOUND
        assert "Pull failed" in caplog.text

    def test_pull_init_that_registers_nothing_fails(self, cli, cli_rg, caplog):
        with patch.object(type(cli_rg.genome), "init_from_remote", return_value=None):
            assert cli("pull", "-g", "nosuchgenome", "--init") == EXIT_NOT_FOUND
        assert "Init failed" in caplog.text

    def test_pull_init_that_registers_something_succeeds(self, cli, cli_rg):
        with patch.object(TransferManager, "init_genomes", return_value=["a" * 32]):
            assert cli("pull", "-g", "rCRSd", "--init") == 0

    def test_pull_init_takes_a_genome_digest(self, cli, cli_rg):
        digest = "d" * 32
        with patch.object(TransferManager, "init_genomes", return_value=[digest]) as init:
            assert cli("pull", "--genome-digest", digest, "--init") == 0
        assert init.call_args.kwargs == {"genomes": [GenomeDigest(digest)], "all_genomes": False}

    @pytest.mark.parametrize("mode", [["--all"], ["--asset", "fasta"]])
    def test_bulk_pull_that_pulls_something_succeeds(self, cli, cli_rg, mode):
        result = BulkPullResult(pulled=[MagicMock()])
        with patch.object(TransferManager, "pull_genomes", return_value=result):
            assert cli("pull", "-g", "rCRSd", *mode) == 0

    @pytest.mark.parametrize("mode", [["--all"], ["--asset", "fasta"]])
    def test_bulk_pull_that_skips_everything_reports_skipped(self, cli, cli_rg, mode, caplog):
        result = BulkPullResult(skipped=[{"asset_group_name": "fasta"}])
        with patch.object(TransferManager, "pull_genomes", return_value=result):
            assert cli("pull", "-g", "rCRSd", *mode) == 0
        assert "Pull skipped" in caplog.text
        assert "Pull failed" not in caplog.text


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
        cli_rg.build.run(
            recipe_name="fasta",
            genome_alias="rCRSd",
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

    def test_id_by_genome_digest(self, cli, cli_rg, capsys):
        """`--genome-digest` names the genome: alone for the genome, with an
        asset path for the asset."""
        digest = cli_rg.alias.resolve("rCRSd")
        assert cli("id", "rCRSd/fasta") == 0
        asset_digest = capsys.readouterr().out.strip()
        assert cli("id", "--genome-digest", digest) == 0
        assert capsys.readouterr().out.strip() == digest
        assert cli("id", "fasta", "--genome-digest", digest) == 0
        assert capsys.readouterr().out.strip() == asset_digest

    def test_getseq_by_genome_digest(self, cli, cli_rg, capsys):
        """`getseq --genome-digest <digest>` works like `getseq -g <alias>`;
        `-g <digest>` is an alias lookup and fails."""
        digest = cli_rg.alias.resolve("rCRSd")
        assert cli("getseq", "-g", "rCRSd", "-l", "rCRSd:1-10") == 0
        by_alias = capsys.readouterr().out.strip()
        assert cli("getseq", "--genome-digest", digest, "-l", "rCRSd:1-10") == 0
        assert capsys.readouterr().out.strip() == by_alias
        assert cli("getseq", "-g", digest, "-l", "rCRSd:1-10") != 0

    def test_build_takes_an_alias_only(self, cli, cli_rg):
        """`build <digest>/fasta` is an alias lookup: the build folder is named
        after the alias, so a digest there is an unknown alias."""
        digest = cli_rg.alias.resolve("rCRSd")
        assert cli("build", f"{digest}/fasta:bydigest") != 0
        assert not (cli_rg.genome_folder / "alias" / "rCRSd" / "fasta" / "bydigest").exists()

    def test_seek_by_genome_digest(self, cli, cli_rg, capsys):
        """`seek fasta --genome-digest <d>` finds what `seek rCRSd/fasta` finds."""
        digest = cli_rg.alias.resolve("rCRSd")
        assert cli("seek", "rCRSd/fasta") == 0
        by_alias = capsys.readouterr().out.strip()
        assert cli("seek", "fasta", "--genome-digest", digest) == 0
        assert capsys.readouterr().out.strip() == by_alias

    def test_genome_in_path_and_flag_is_an_error(self, cli, cli_rg):
        """The genome is named once: the path's alias or the flag, not both."""
        digest = cli_rg.alias.resolve("rCRSd")
        assert cli("seek", "rCRSd/fasta", "--genome-digest", digest) == 3
        assert cli("remove", "rCRSd/fasta", "-f", "--genome-digest", digest) == 3

    def test_malformed_genome_digest_is_invalid_input(self, cli):
        assert cli("seek", "fasta", "--genome-digest", "rCRSd") == 3

    def test_a_genome_with_no_alias_is_reachable_by_digest(self, cli, cli_rg, capsys):
        """With its only alias gone, the genome is still reachable by digest."""
        digest = cli_rg.alias.resolve("rCRSd")
        cli_rg.alias.remove(GenomeAlias("rCRSd"))
        assert cli("seek", "rCRSd/fasta") == 2
        assert cli("seek", "fasta", "--genome-digest", digest) == 0
        assert f"/data/{digest}/fasta/" in capsys.readouterr().out

    def test_alias_spelled_like_a_digest(self, cli, cli_rg, capsys):
        """An alias spelled like another genome's digest is only ever an alias
        in a registry path, and that digest only ever a digest via the flag."""
        rcrsd = cli_rg.alias.resolve("rCRSd")
        other = fake_digest("other genome")
        cli_rg.genome.add(other, "a genome with no assets", [])
        cli_rg.alias.add(GenomeAlias(other), rcrsd)

        # Through the path: the alias, which names rCRSd, which has a fasta asset.
        capsys.readouterr()
        assert cli("seek", f"{other}/fasta", "--abs") == 0
        assert f"/data/{rcrsd}/fasta/" in capsys.readouterr().out
        # Through the flag: the digest, whose genome has no fasta asset.
        assert cli("seek", "fasta", "--genome-digest", other) == 2

    def test_list_by_genome_digest(self, cli, cli_rg, capsys):
        digest = cli_rg.alias.resolve("rCRSd")
        assert cli("list", "--genome-digest", digest) == 0
        out = capsys.readouterr().out
        assert "fasta" in out and "default" in out
        assert cli("list", "-g", digest) != 0

    def test_genome_remove_by_genome_digest(self, cli, cli_rg):
        digest = cli_rg.alias.resolve("rCRSd")
        assert cli("genome", "remove", "--genome", digest, "-f") != 0
        assert cli_rg.genome.exists(digest)
        assert cli("genome", "remove", "--genome-digest", digest, "-f") == 0
        assert not cli_rg.genome.exists(digest)

    def test_pull_passes_the_typed_genome(self, cli, cli_rg):
        """pull hands the library a GenomeAlias from a path, a GenomeDigest from the flag."""
        digest = fake_digest("remote genome")
        with (
            patch.object(TransferManager, "pull", return_value=MagicMock()) as pull,
            patch("refgenie.cli.commands.pull.genome_digest_from", return_value=digest),
            patch("refgenie.managers.asset.manager.AssetManager.seek", return_value="/x"),
        ):
            assert cli("pull", "hg38/fasta") == 0
            assert cli("pull", "fasta", "--genome-digest", digest) == 0
        by_alias, by_digest = (c.kwargs["genome"] for c in pull.call_args_list)
        assert type(by_alias) is GenomeAlias and by_alias == "hg38"
        assert type(by_digest) is GenomeDigest and by_digest == digest

    def test_build_requirements(self, cli, cli_rg_with_recipe, capsys):
        assert cli("build", "rCRSd/fasta", "--requirements") == 0
        assert "fasta" in capsys.readouterr().out.lower()

    def test_alias_set_no_genome_errors(self, cli):
        """alias set with a nonexistent digest exits non-zero, unless forced."""
        missing = "n" * 32
        assert cli("alias", "set", "-a", "my_alias", "-d", missing) != 0
        assert cli("alias", "set", "-a", "my_alias", "-d", missing, "--force") == 0
        # A malformed digest is invalid input, forced or not.
        assert cli("alias", "set", "-a", "other", "-d", "nonexistent_digest", "--force") == 3

    def test_remove_aliases_drops_alias_when_last_asset_goes(self, cli, cli_engine):
        """`remove --aliases` leaves no alias behind once nothing is left."""
        assert cli("remove", "rCRSd/fasta", "-f", "--aliases") == 0
        fresh = Refgenie(database_engine=cli_engine, suppress_migrations=True)
        assert [a.name for a in fresh.alias.list_all()] == []


# ---------------------------------------------------------------------------
# serve / dash / generate dispatch (mocked servers)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "argv, module, runner, expected",
    [
        (["serve"], "refgenie.server.main", "run_server", {"port": 8000, "reload": False}),
        (["serve", "--port", "9999"], "refgenie.server.main", "run_server", {"port": 9999}),
        (["serve", "--reload"], "refgenie.server.main", "run_server", {"reload": True}),
        (["dash"], "refgenie.server.main", "run_local", {"port": 8080}),
        (["dash", "--port", "7777"], "refgenie.server.main", "run_local", {"port": 7777}),
        (
            ["generate", "snakefile", "-o", "/tmp/Snakefile"],
            "refgenie.integrations.snakemake.generate",
            "populate_snakefile_template",
            {"snakefile_output_path": Path("/tmp/Snakefile")},
        ),
    ],
    ids=["serve", "serve-port", "serve-reload", "dash", "dash-port", "generate-snakefile"],
)
def test_command_dispatches_to_its_runner(cli, cli_rg, argv, module, runner, expected):
    """serve, dash and generate reach their runner with the parsed arguments."""
    mock_runner = MagicMock()
    with (
        patch.dict("sys.modules", {module: MagicMock(**{runner: mock_runner})}),
        patch("webbrowser.open_new_tab") as open_tab,
    ):
        assert cli(*argv) == 0
    kwargs = mock_runner.call_args.kwargs
    for name, value in expected.items():
        assert kwargs[name] == value
    if argv[0] == "dash":
        open_tab.assert_called_once_with(f"http://127.0.0.1:{expected['port']}")


# ---------------------------------------------------------------------------
# `refgenie id` and the remote-digest lookup behind it
# ---------------------------------------------------------------------------


@pytest.fixture
def cli_genome(cli_rg, fixtures_path):
    """``cli_rg`` with the rCRSd genome initialized (nothing built)."""
    cli_rg.genome.initialize_genome(
        fasta_file_path=fixtures_path / "rCRSd.fa", alias_names=["rCRSd"], description="rCRSd"
    )
    return cli_rg


@pytest.mark.component
class TestIdFlags:
    def test_verbose_output_contains_metadata(self, cli, cli_genome, capsys):
        assert cli("id", "rCRSd", "--verbose") == 0
        out = capsys.readouterr().out
        for key in ("digest:", "sequences:", "total_length:", "source:"):
            assert key in out

    def test_validate_succeeds_for_valid_genome(self, cli, cli_genome, capsys):
        assert cli("id", "rCRSd", "--validate-store") == 0
        assert capsys.readouterr().out.strip()

    def test_info_shows_aliases_and_metadata(self, cli, cli_genome, capsys):
        digest = cli_genome.alias.resolve("rCRSd")
        assert cli("id", digest, "--info") == 0
        out = capsys.readouterr().out
        for key in ("digest:", "aliases:", "rCRSd", "sequences:", "total_length:"):
            assert key in out

    def test_info_unknown_digest_errors(self, cli, cli_genome):
        assert cli("id", "nonexistent_digest_0000000000000", "--info") != 0


@pytest.mark.component
class TestIdRemote:
    def test_remote_with_unknown_alias_suggests_genome_digest(self, cli, cli_genome, caplog):
        """A positional name is an alias, even when it is shaped like a digest."""
        with patch.object(ServerManager, "locate") as locate:
            assert cli("id", "a" * 32, "--remote") == EXIT_INVALID_INPUT
        assert "--genome-digest" in caplog.text
        locate.assert_not_called()

    def test_remote_with_digest_asks_the_servers(self, cli, cli_genome, capsys):
        digest = "a" * 32
        with patch.object(ServerManager, "locate", return_value=("http://s", digest, {})) as loc:
            assert cli("id", "--genome-digest", digest, "--remote") == 0
        loc.assert_called_once_with(digest)
        assert isinstance(loc.call_args.args[0], GenomeDigest)
        assert digest in capsys.readouterr().out

    def test_remote_with_digest_not_found_errors(self, cli, cli_genome):
        with patch.object(ServerManager, "locate", return_value=None):
            assert cli("id", "--genome-digest", "b" * 32, "--remote") == EXIT_NOT_FOUND

    def test_malformed_digest_is_invalid_input(self, cli, cli_genome):
        assert cli("id", "--genome-digest", "not-a-digest", "--remote") == EXIT_INVALID_INPUT

    def test_known_digest_is_found_locally(self, cli, cli_genome, capsys):
        digest = cli_genome.alias.resolve("rCRSd")
        with patch.object(ServerManager, "locate") as locate:
            assert cli("id", "--genome-digest", digest, "--remote") == 0
        assert capsys.readouterr().out.strip() == digest
        locate.assert_not_called()

    def test_unknown_digest_without_remote_errors(self, cli, cli_genome):
        assert cli("id", "--genome-digest", "c" * 32) == EXIT_NOT_FOUND


REMOTE_DIGEST = GenomeDigest("r" * 32)
SERVER = "http://server.test"


# Two genomes with the same two sequences under different names.
_SIDES = (
    {"names": ["chr1", "chr2"], "lengths": [10, 20], "sequences": ["SQ.a", "SQ.b"]},
    {"names": ["1", "2"], "lengths": [10, 20], "sequences": ["SQ.a", "SQ.b"]},
)


class TestRemoteGenomeByName:
    """Commands that reach a genome only a server knows, through ``locate``."""

    def test_listr_finds_an_alias_on_the_servers(self, cli, cli_rg):
        with (
            patch.object(ServerManager, "locate", return_value=(SERVER, REMOTE_DIGEST, {})) as loc,
            patch.object(ServerManager, "assets_table", return_value=[]) as table,
        ):
            assert cli("listr", "-g", "hg38") == 0
        assert loc.call_args.args == (GenomeAlias("hg38"), None)
        assert table.call_args.kwargs["genome_digests"] == [REMOTE_DIGEST]

    def test_listr_of_a_genome_no_one_knows_is_not_found(self, cli, cli_rg):
        with (
            patch.object(ServerManager, "locate", return_value=None),
            patch.object(ServerManager, "assets_table", return_value=[]) as table,
        ):
            assert cli("listr", "-g", "hg38") == EXIT_NOT_FOUND
        table.assert_not_called()

    def test_getseq_registers_the_genome_first(self, cli, cli_rg, capsys):
        from refgenie.managers.sequence import SequenceManager

        with (
            patch.object(
                type(cli_rg.genome), "init_from_remote", return_value=REMOTE_DIGEST
            ) as init,
            patch.object(SequenceManager, "get", return_value="ACGT") as get,
        ):
            assert cli("getseq", "-g", "hg38", "-l", "chr1:1-4") == 0
        assert init.call_args.args == (GenomeAlias("hg38"),)
        assert init.call_args.kwargs == {"require_collection": True}
        assert get.call_args.args == (REMOTE_DIGEST, "chr1", 0, 4)
        assert "ACGT" in capsys.readouterr().out

    def test_getseq_of_a_genome_no_one_can_provide_is_not_found(self, cli, cli_rg):
        with patch.object(type(cli_rg.genome), "init_from_remote", return_value=None):
            assert cli("getseq", "-g", "hg38", "-l", "chr1:1-4") == EXIT_NOT_FOUND

    # Local genomes are named by digest: a bare genome's alias lives in this
    # process's alias cache only, so the CLI's own Refgenie would not see it.

    def test_compare_of_two_local_genomes_asks_no_server(self, cli, cli_rg, capsys):
        cli_rg.genome.add("a" * 32, "a", [])
        cli_rg.genome.add("b" * 32, "b", [])
        with (
            patch.object(type(cli_rg.genome), "compare_sides", return_value=_SIDES) as sides,
            patch.object(ServerManager, "client") as client,
        ):
            assert cli("compare", "--genome-digest", f"{'a' * 32},{'b' * 32}") == 0
        assert sides.call_args.args == ("a" * 32, "b" * 32)
        client.assert_not_called()
        out = capsys.readouterr().out
        assert "Same 2 sequences in the same order; all 2 names differ" in out
        assert "chr1 -> 1" in out

    def test_compare_with_a_remote_genome_asks_its_server(self, cli, cli_rg, capsys):
        cli_rg.genome.add("a" * 32, "a", [])
        from refgenie.utils.seqcol import compare_level2

        client = MagicMock()
        comparison = compare_level2(*_SIDES, REMOTE_DIGEST, "a" * 32)
        client.genome.compare.return_value = {
            **comparison,
            "summary": "worded by the server",
            "name_map": {"chr1": "1", "chr2": "2"},
        }
        with (
            patch.object(ServerManager, "locate", return_value=(SERVER, REMOTE_DIGEST, {})),
            patch.object(ServerManager, "client", return_value=client) as get_client,
        ):
            assert cli("compare", "hg38", "--genome-digest", "a" * 32, "--json") == 0
        get_client.assert_called_once_with(SERVER)
        client.genome.compare.assert_called_once_with(REMOTE_DIGEST, "a" * 32)
        import json

        payload = json.loads(capsys.readouterr().out)
        # Worded here, with the names the user typed.
        assert payload["summary"].startswith("Same 2 sequences in the same order")
        assert payload["name_map"] == {"chr1": "1", "chr2": "2"}
        assert payload["comparison"]["digests"]["a"] == REMOTE_DIGEST

    def test_compare_the_server_cannot_answer_suggests_registering(self, cli, cli_rg, caplog):
        cli_rg.genome.add("a" * 32, "a", [])
        client = MagicMock()
        client.genome.compare.side_effect = ValueError("404")
        with (
            patch.object(ServerManager, "locate", return_value=(SERVER, REMOTE_DIGEST, {})),
            patch.object(ServerManager, "client", return_value=client),
        ):
            assert cli("compare", "hg38", "--genome-digest", "a" * 32) == EXIT_NOT_FOUND
        assert "Register 'hg38' here first (refgenie pull --init" in caplog.text

    def test_genome_init_from_a_server_reads_the_located_store(self, cli, cli_rg):
        import refgenie.cli.commands.genome as handlers_mod

        attrs = {"refget_store_url": "https://store.test/igenomes/"}
        source = MagicMock()
        with (
            patch.object(
                ServerManager, "locate", return_value=(SERVER, REMOTE_DIGEST, attrs)
            ) as loc,
            patch.object(handlers_mod, "collection_source", return_value=source) as pick,
            patch.object(
                BuildManager, "initialize_and_build", return_value=(REMOTE_DIGEST, True)
            ) as init,
        ):
            code = cli(
                "genome",
                "init",
                "--server",
                SERVER,
                "--digest",
                REMOTE_DIGEST,
                "-n",
                "x",
                "--no-build",
            )
        assert code == 0
        assert loc.call_args.args == (GenomeDigest(REMOTE_DIGEST), [SERVER])
        pick.assert_called_once_with(SERVER, attrs)
        assert init.call_args.kwargs["source"] is source
        assert init.call_args.kwargs["digest"] == REMOTE_DIGEST

    def test_genome_init_from_a_server_that_lacks_the_genome(self, cli, cli_rg):
        with patch.object(ServerManager, "locate", return_value=None):
            code = cli("genome", "init", "--server", SERVER, "-n", "x", "--no-build")
        assert code == EXIT_NOT_FOUND


# ---------------------------------------------------------------------------
# `refgenie-build-fasta`: FASTA, .fai and .chrom.sizes out of a RefgetStore
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
    def test_writes_fasta_index_and_chrom_sizes(self, store_with_genome, tmp_path):
        from gtars.refget import RefgetStore
        from refgenie.cli.build_fasta import build_fasta_main

        store_path, digest = store_with_genome
        output = tmp_path / "output"
        build_fasta_main([str(store_path), digest, str(output)])

        fa_content = (output / f"{digest}.fa").read_text()
        assert fa_content.startswith(">") and len(fa_content) > 100

        fai_lines = (output / f"{digest}.fa.fai").read_text().strip().split("\n")
        assert fai_lines and all(len(line.split("\t")) == 5 for line in fai_lines)

        level2 = RefgetStore.open_local(store_path).get_collection_level2(digest)
        lines = (output / f"{digest}.chrom.sizes").read_text().strip().split("\n")
        assert [line.split("\t") for line in lines] == [
            [name, str(length)] for name, length in zip(level2["names"], level2["lengths"])
        ]

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


# ===========================================================================
# Regression guard: every CLI model field is consumed by its handler
# ===========================================================================

# Fields consumed indirectly, not via ``cmd.<field>`` in the handler body.
# Keep this short: growth here means a flag went unwired again.
_INDIRECT = {
    # LargeArchiveFlags reads these through ``self.<field>`` in
    # resolve_force_large() and its validators (refgenie/cli/commands/pull.py).
    "skip_large",
    "pull_large",
    "batch",
}


def _model_handler_pairs():
    """Yield (model, handler) for every leaf command model in the CLI."""
    pairs = []
    for model, handler in get_dispatch().items():
        group_dispatch = None
        for name, value in vars(inspect.getmodule(handler)).items():
            if name.endswith("_DISPATCH") and isinstance(value, dict):
                group_dispatch = value
        if model.__name__.endswith("Parser"):
            if group_dispatch:
                pairs.extend(group_dispatch.items())
        else:
            pairs.append((model, handler))
    # asset.handle_asset_group forwards its only leaf to listing.handle_list.
    pairs.append((cli_asset.AssetListNestedModel, cli_listing.handle_list))
    return pairs


_PAIRS = _model_handler_pairs()


@pytest.mark.parametrize(
    "model,handler",
    _PAIRS,
    ids=[f"{m.__name__}-{h.__name__}" for m, h in _PAIRS],
)
def test_every_model_field_is_consumed(model, handler):
    """A field with no handler reader renders a flag that is a silent no-op.
    The rule is wire-or-delete: if this fails, consume the field or remove it --
    do not extend ``_INDIRECT``."""
    source = inspect.getsource(handler)
    # A handler may delegate ``cmd`` to a shared helper or to a helper in its
    # own module (which may delegate again); count what those read as consumed.
    candidates = {
        name: value
        for module in (cli_helpers, inspect.getmodule(handler))
        for name, value in vars(module).items()
        if inspect.isfunction(value) and value is not handler
    }
    followed = set()
    while new := [n for n in candidates if n not in followed and re.search(rf"\b{n}\(", source)]:
        for name in new:
            followed.add(name)
            source += inspect.getsource(candidates[name])
    unwired = [
        field
        for field in model.model_fields
        if field not in _INDIRECT
        and not re.search(rf"cmd\.{field}\b", source)
        and not re.search(rf"""getattr\(\s*cmd\s*,\s*["']{field}["']""", source)
    ]
    assert not unwired, (
        f"{model.__name__} fields {unwired} are never read by "
        f"{handler.__module__}.{handler.__name__}. Wire the flag or delete the field."
    )


class TestPluginsCommand:
    """`refgenie plugins`: list installed plugins and manage their settings."""

    def test_set_then_list_shows_the_setting(self, cli, cli_rg, capsys):
        assert cli("plugins", "set", "nfcore", "config_path=/x/nf.config") == 0
        assert cli_rg.plugins.settings("nfcore") == {"config_path": "/x/nf.config"}
        capsys.readouterr()
        assert cli("plugins") == 0
        out = capsys.readouterr().out
        assert "config_path" in out
        assert "/x/nf.config" in out
        assert "(not installed)" in out

    def test_unset_removes_the_setting(self, cli, cli_rg):
        cli_rg.plugins.set("nfcore", config_path="/x/nf.config", other="y")
        assert cli("plugins", "unset", "nfcore", "config_path") == 0
        assert cli_rg.plugins.settings("nfcore") == {"other": "y"}
        assert cli("plugins", "unset", "nfcore") == 0
        assert cli_rg.plugins.all_settings() == {}

    def test_set_without_equals_fails(self, cli, cli_rg):
        assert cli("plugins", "set", "nfcore", "bad") != 0
        assert cli_rg.plugins.settings("nfcore") == {}

    def test_lists_an_installed_plugin(self, cli, cli_rg, capsys, install_plugins):
        install_plugins({"post_update": [("fakeplugin", "tests.plugins.fake_plugin:record")]})
        assert cli("plugins", "list") == 0
        out = capsys.readouterr().out
        assert "post_update" in out
        assert "fakeplugin" in out
        assert "ok" in out

    def test_no_plugins_installed(self, cli, cli_rg, capsys, install_plugins):
        install_plugins({})
        assert cli("plugins") == 0
        assert "No refgenie plugins installed" in capsys.readouterr().out

    def test_one_command_fires_one_post_update(self, cli, cli_rg, install_plugins):
        """A command that runs several operations (like `genome init --fasta`,
        which initializes then builds) still fires `post_update` once."""
        from refgenie.plugins.hooks import Change
        from tests.plugins import fake_plugin

        install_plugins({"post_update": [("fakeplugin", "tests.plugins.fake_plugin:record")]})

        def two_operations(cmd, rg):
            rg._events.record(Change("genome_added", genome="a"))
            rg._events.record(Change("asset_added", genome="a"))

        with patch("refgenie.cli.dispatch.get_dispatch") as get:
            get.return_value.get.return_value = two_operations
            assert cli("list") == 0
        updates = [event for hook, event in fake_plugin.CALLS if hook == "post_update"]
        assert len(updates) == 1
        assert [c.action for c in updates[0].changes] == ["genome_added", "asset_added"]


class TestQuietScriptableOutput:
    """getseq, --json listings, YAML config and the log level."""

    @pytest.fixture(autouse=True)
    def setup_asset(self, cli_rg, fixtures_path):
        register_fasta(cli_rg, fixtures_path)
        cli_rg.genome.initialize_genome(
            fasta_file_path=fixtures_path / "rCRSd.fa",
            alias_names=["rCRSd"],
            description="rCRSd mitochondrial reference",
        )
        cli_rg.build.run(
            recipe_name="fasta",
            genome_alias="rCRSd",
            asset_group_name="fasta",
            asset_name="default",
        )

    def test_getseq_plain_is_one_line(self, cli, capsys):
        assert cli("getseq", "-g", "rCRSd", "-l", "rCRSd") == 0
        out = capsys.readouterr().out
        assert out.count("\n") == 1 and len(out.strip()) == 33138

    def test_getseq_fasta(self, cli, capsys):
        assert cli("getseq", "-g", "rCRSd", "-l", "rCRSd:1-10", "--fasta") == 0
        header, bases = capsys.readouterr().out.splitlines()
        assert header == ">rCRSd:1-10" and len(bases) == 10
        assert cli("getseq", "-g", "rCRSd", "-l", "rCRSd", "--fasta") == 0
        header, *lines = capsys.readouterr().out.splitlines()
        assert header == ">rCRSd"
        assert all(len(x) == 60 for x in lines[:-1]) and sum(map(len, lines)) == 33138

    def test_getseq_zero_based(self, cli, capsys):
        """The same 10 bases, as a 1-based text region or with --zero-based;
        stderr says how the region was read."""
        assert cli("getseq", "-g", "rCRSd", "-l", "rCRSd:1-10") == 0
        one = capsys.readouterr()
        assert cli("getseq", "-g", "rCRSd", "-l", "rCRSd:0-10", "--zero-based") == 0
        assert capsys.readouterr().out == one.out
        assert len(one.out.strip()) == 10
        assert "rCRSd:1-10 (1-based, inclusive): 10 bases" in one.err

    def test_getseq_unknown_name(self, cli):
        assert cli("getseq", "-g", "rCRSd", "-l", "chr999") == EXIT_NOT_FOUND

    def test_list_json(self, cli, capsys):
        import json

        assert cli("list", "--json") == 0
        data = json.loads(capsys.readouterr().out)
        assert data[0]["source"] == "local"
        genome = data[0]["genomes"][0]
        assert genome["aliases"] == ["rCRSd"] and "fasta:default" in genome["assets"]

    def test_listr_json(self, cli, capsys):
        import json

        listing = ({"http://s": {"d1": ["fasta:default"]}}, {"http://s": {"d1": ["a"]}})
        with patch.object(ServerManager, "list_assets", return_value=listing):
            assert cli("listr", "--json") == 0
        data = json.loads(capsys.readouterr().out)
        assert data == [
            {
                "source": "http://s",
                "genomes": [{"digest": "d1", "aliases": ["a"], "assets": ["fasta:default"]}],
            }
        ]

    def test_config_get_is_yaml(self, cli, capsys):
        import yaml

        assert cli("config", "get") == 0
        data = yaml.safe_load(capsys.readouterr().out)
        assert set(data) == {"configuration", "environment"}

    def test_info_logs_hidden_but_reach_logger_handlers(self, capsys):
        import logging

        from refgenie.logger import logger

        seen = []

        class _H(logging.Handler):
            def emit(self, record):
                seen.append(record.getMessage())

        h = _H(logging.INFO)
        logger.addHandler(h)
        try:
            logger.info("quiet please")
        finally:
            logger.removeHandler(h)
        assert seen == ["quiet please"]
        assert "quiet please" not in capsys.readouterr().err
