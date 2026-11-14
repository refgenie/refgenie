"""Parse argv the way the refgenie CLI does, without running a command."""


def parse(argv):
    """Parse argv as ``main`` does."""
    from refgenie.cli.main import _make_cli_source
    from refgenie.cli.parser import TopLevelParser

    source = _make_cli_source(TopLevelParser, argv)
    return TopLevelParser(_cli_settings_source=source)


def subcmd(argv):
    """Parse argv and return the active leaf subcommand model."""
    from pydantic_settings import get_subcommand

    return get_subcommand(parse(argv), is_required=True)
