"""Tests for refgenie.argparser module."""

from __future__ import annotations

import pytest

from refgenie.argparser import build_argparser
from refgenie.const import (
    ALIAS_CMD,
    BUILD_CMD,
    COMPARE_CMD,
    GET_ASSET_CMD,
    GET_REMOTE_ASSET_CMD,
    GETSEQ_CMD,
    ID_CMD,
    INIT_CMD,
    INSERT_CMD,
    LIST_LOCAL_CMD,
    LIST_REMOTE_CMD,
    POPULATE_CMD,
    POPULATE_REMOTE_CMD,
    PULL_CMD,
    REMOVE_CMD,
    SUBSCRIBE_CMD,
    TAG_CMD,
    UNSUBSCRIBE_CMD,
    UPGRADE_CMD,
)


@pytest.fixture
def parser():
    return build_argparser()


class TestBuildArgparser:
    """Tests for argparser construction and parsing."""

    def test_returns_parser(self, parser):
        import argparse

        assert isinstance(parser, argparse.ArgumentParser)

    def test_init_command(self, parser):
        args = parser.parse_args(["init", "-c", "/path/to/config.yaml"])
        assert args.command == INIT_CMD
        assert args.genome_config == "/path/to/config.yaml"

    def test_list_command(self, parser):
        args = parser.parse_args(["list", "-c", "g.yaml"])
        assert args.command == LIST_LOCAL_CMD

    def test_build_command_with_asset(self, parser):
        args = parser.parse_args(["build", "-c", "g.yaml", "hg38/fasta"])
        assert args.command == BUILD_CMD
        assert args.asset_registry_paths == ["hg38/fasta"]

    def test_build_reduce_flag(self, parser):
        args = parser.parse_args(["build", "-c", "g.yaml", "--reduce"])
        assert args.command == BUILD_CMD
        assert args.reduce is True

    def test_seek_command(self, parser):
        args = parser.parse_args(["seek", "-c", "g.yaml", "hg38/fasta"])
        assert args.command == GET_ASSET_CMD
        assert args.asset_registry_paths == ["hg38/fasta"]

    def test_pull_command(self, parser):
        args = parser.parse_args(["pull", "-c", "g.yaml", "hg38/fasta"])
        assert args.command == PULL_CMD

    def test_pull_batch_mode(self, parser):
        args = parser.parse_args(["pull", "-c", "g.yaml", "-b", "hg38/fasta"])
        assert args.batch is True

    def test_id_command(self, parser):
        args = parser.parse_args(["id", "-c", "g.yaml", "hg38/fasta:default"])
        assert args.command == ID_CMD

    def test_subscribe_command(self, parser):
        args = parser.parse_args(
            ["subscribe", "-c", "g.yaml", "-s", "http://server.com"]
        )
        assert args.command == SUBSCRIBE_CMD
        assert args.genome_server == ["http://server.com"]

    def test_compare_command(self, parser):
        args = parser.parse_args(["compare", "-c", "g.yaml", "hg38", "mm10"])
        assert args.command == COMPARE_CMD
        assert args.genome1 == ["hg38"]
        assert args.genome2 == ["mm10"]

    def test_all_subcommands_registered(self, parser):
        """Verify all expected commands are registered as subparsers."""
        subparsers_actions = [
            a for a in parser._subparsers._actions if hasattr(a, "_parser_class")
        ]
        assert len(subparsers_actions) == 1
        registered = set(subparsers_actions[0].choices.keys())
        expected = {
            INIT_CMD,
            LIST_LOCAL_CMD,
            LIST_REMOTE_CMD,
            PULL_CMD,
            BUILD_CMD,
            GET_ASSET_CMD,
            GET_REMOTE_ASSET_CMD,
            INSERT_CMD,
            REMOVE_CMD,
            GETSEQ_CMD,
            TAG_CMD,
            ID_CMD,
            SUBSCRIBE_CMD,
            UNSUBSCRIBE_CMD,
            ALIAS_CMD,
            COMPARE_CMD,
            UPGRADE_CMD,
            POPULATE_CMD,
            POPULATE_REMOTE_CMD,
        }
        assert expected.issubset(registered), (
            f"Missing commands: {expected - registered}"
        )

    def test_no_command_gives_none(self, parser):
        args, _ = parser.parse_known_args([])
        assert args.command is None
