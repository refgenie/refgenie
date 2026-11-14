"""What `refgenie seek` tells a user when the asset is not here."""

import pytest

from refgenie.cli.commands.helpers import _lookup_one
from refgenie.cli.errors import EXIT_NOT_FOUND
from refgenie.core.root import Refgenie
from refgenie.exceptions import MissingAssetGroupError, MissingSeekKeyError

pytestmark = pytest.mark.unit


class _Parser:
    parse_asset_registry_path = staticmethod(Refgenie.parse_asset_registry_path)


def _raise(exc):
    def lookup(parsed):
        raise exc

    return lookup


def test_a_missing_asset_names_the_pull_command(capsys):
    code = _lookup_one(
        _Parser(),
        "R64-1-1/bwa_index.dir",
        _raise(MissingAssetGroupError("d", "bwa_index")),
        pull_hint=True,
    )
    assert code == EXIT_NOT_FOUND
    assert "refgenie pull R64-1-1/bwa_index" in capsys.readouterr().err


def test_a_missing_seek_key_gets_no_pull_hint(capsys):
    code = _lookup_one(
        _Parser(),
        "R64-1-1/fasta.nokey",
        _raise(MissingSeekKeyError("d", "fasta", "default", "nokey")),
        pull_hint=True,
    )
    assert code == EXIT_NOT_FOUND
    assert "refgenie pull" not in capsys.readouterr().err
