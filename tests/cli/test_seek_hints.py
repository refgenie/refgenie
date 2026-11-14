"""What `refgenie seek` tells a user when the asset is not here."""

from types import SimpleNamespace

import pytest

from refgenie.cli.commands.helpers import _lookup_one
from refgenie.cli.errors import EXIT_NOT_FOUND
from refgenie.core.root import Refgenie
from refgenie.exceptions import MissingAssetGroupError, MissingSeekKeyError

pytestmark = pytest.mark.unit

DIGEST = "d" * 32


def _refgenie(pulled: list[str]):
    """A stand-in with hg38 pulled here, holding ``pulled`` (``group.key:asset``)."""
    return SimpleNamespace(
        parse_asset_registry_path=Refgenie.parse_asset_registry_path,
        alias=SimpleNamespace(resolve=lambda name: DIGEST),
        asset=SimpleNamespace(list_all=lambda **kw: ({DIGEST: pulled}, {})),
    )


def _raise(exc):
    def lookup(parsed):
        raise exc

    return lookup


FASTA_INDEX = ["fasta_index.fai:default", "fasta_index.chrom_sizes:default"]


def test_a_seek_key_named_like_the_group_is_suggested(capsys):
    code = _lookup_one(
        _refgenie(FASTA_INDEX),
        "hg38/chrom_sizes",
        _raise(MissingAssetGroupError(DIGEST, "chrom_sizes")),
        pull_hint=True,
    )
    err = capsys.readouterr().err
    assert code == EXIT_NOT_FOUND
    assert "Did you mean: hg38/fasta_index.chrom_sizes?" in err
    assert "refgenie pull" not in err


def test_an_asset_not_here_names_the_pull_command(capsys):
    code = _lookup_one(
        _refgenie(FASTA_INDEX),
        "hg38/bwa_index.dir",
        _raise(MissingAssetGroupError(DIGEST, "bwa_index")),
        pull_hint=True,
    )
    err = capsys.readouterr().err
    assert code == EXIT_NOT_FOUND
    assert "refgenie pull hg38/bwa_index" in err
    assert "Did you mean" not in err


def test_a_misspelled_seek_key_gets_a_suggestion_and_no_pull(capsys):
    code = _lookup_one(
        _refgenie(FASTA_INDEX),
        "hg38/fasta_index.chromsizes",
        _raise(MissingSeekKeyError(DIGEST, "fasta_index", "default", "chromsizes")),
        pull_hint=True,
    )
    err = capsys.readouterr().err
    assert code == EXIT_NOT_FOUND
    assert "hg38/fasta_index.chrom_sizes" in err
    assert "refgenie pull" not in err
