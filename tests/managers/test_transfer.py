"""
Unit tests for TransferManager's bulk-pull skeleton (``rgc.transfer``).

Every dependency is a MagicMock: no database, no server. These pin the shared
collect -> confirm -> pull_multiple path that pull_genomes and mirror use, and
the alias-list handling in init_genomes.
"""

import logging
from unittest.mock import MagicMock

import pytest

from refgenie.exceptions import MissingAliasError
from refgenie.managers.transfer import TransferManager
from refgenie.models import GenomeAlias, GenomeDigest
from refgenie.plugins.events import EventSink

pytestmark = pytest.mark.unit

D1 = GenomeDigest("a" * 32)
D2 = GenomeDigest("b" * 32)


def _asset(group: str, digest: str = D1, size: int = 10) -> dict:
    return {"asset_group_name": group, "genome_digest": digest, "archive_size": size}


@pytest.fixture
def deps():
    servers = MagicMock()
    alias = MagicMock()
    genome = MagicMock()
    puller = MagicMock()
    puller.pull_multiple.side_effect = lambda asset_list, **kwargs: list(asset_list)
    return servers, alias, genome, puller


@pytest.fixture
def transfer(deps):
    servers, alias, genome, puller = deps
    return TransferManager(
        servers=servers,
        alias_manager=alias,
        genome_manager=genome,
        puller=puller,
        events=EventSink(),
    )


def test_pull_genomes_filters_by_asset_group(transfer, deps):
    servers, _, _, puller = deps
    servers.list_assets_for_genome.return_value = [_asset("fasta"), _asset("bowtie2_index")]

    pulled = transfer.pull_genomes([D1], asset_group_name="fasta", force=True)

    assert pulled == [_asset("fasta")]
    assert puller.pull_multiple.call_args.kwargs["asset_list"] == [_asset("fasta")]


def test_pull_genomes_without_group_pulls_everything(transfer, deps):
    servers, _, _, puller = deps
    servers.list_assets_for_genome.return_value = [_asset("fasta"), _asset("bowtie2_index")]

    assert len(transfer.pull_genomes([D1], force=True)) == 2
    puller.pull_multiple.assert_called_once()


def test_pull_genomes_all_genomes_expands_from_servers(transfer, deps):
    servers, alias, _, _ = deps
    servers.list_genomes.return_value = [
        {"aliases": ["hg38"], "genome_digest": D1},
        {"aliases": [], "genome_digest": D2},
    ]
    alias.resolve.return_value = D1
    servers.list_assets_for_genome.side_effect = lambda d: [_asset("fasta", d)]

    pulled = transfer.pull_genomes(all_genomes=True, asset_group_name="fasta", force=True)

    # The aliased genome resolves through the alias manager; the alias-less one
    # is used by digest directly.
    alias.resolve.assert_called_once_with(GenomeAlias("hg38"))
    asked = [call.args[0] for call in servers.list_assets_for_genome.call_args_list]
    assert asked == [D1, D2]
    assert isinstance(asked[1], GenomeDigest)
    assert len(pulled) == 2


def test_refused_confirmation_pulls_nothing(transfer, deps, caplog):
    servers, _, _, puller = deps
    servers.list_assets_for_genome.return_value = [_asset("fasta")]

    with caplog.at_level(logging.INFO, logger="refgenie"):
        assert transfer.pull_genomes([D1], confirm=lambda _msg: False) == []

    puller.pull_multiple.assert_not_called()
    assert "Bulk pull cancelled" in caplog.text


def test_no_assets_warns_and_returns_empty(transfer, deps, caplog):
    servers, _, _, puller = deps
    servers.list_assets_for_genome.return_value = [_asset("bowtie2_index")]

    with caplog.at_level(logging.WARNING, logger="refgenie"):
        assert transfer.pull_genomes([D1], asset_group_name="fasta", force=True) == []

    puller.pull_multiple.assert_not_called()
    assert "No remote assets named 'fasta' found for specified genomes" in caplog.text


def test_init_genomes_all_skips_aliasless_remote_genomes(transfer, deps):
    servers, alias, genome, _ = deps
    servers.list_genomes.return_value = [
        {"aliases": ["hg38"], "genome_digest": D1},
        {"aliases": [], "genome_digest": D2},
    ]
    genome.init_from_remote.return_value = True
    alias.resolve.return_value = D1

    assert transfer.init_genomes(all_genomes=True) == [D1]
    genome.init_from_remote.assert_called_once_with(GenomeAlias("hg38"))


def test_init_genomes_skips_alias_that_does_not_resolve(transfer, deps):
    _, alias, genome, _ = deps
    genome.init_from_remote.return_value = True
    alias.resolve.side_effect = [D1, MissingAliasError("gone")]

    assert transfer.init_genomes([GenomeAlias("hg38"), GenomeAlias("mm10")]) == [D1]


def test_mirror_registers_each_genome_with_description_and_digest_fallback(transfer, deps):
    servers, _, genome, puller = deps
    servers.list_genomes.return_value = [
        {"aliases": ["hg38"], "genome_digest": D1, "description": "human"},
        {"aliases": [], "genome_digest": D2, "description": "no alias"},
    ]
    servers.list_assets_for_genome.side_effect = lambda d: [_asset("fasta", d)]

    pulled = transfer.mirror(force=True)

    registered = [call.kwargs for call in genome.init_from_remote.call_args_list]
    assert registered == [
        {
            "alias_name": GenomeAlias("hg38"),
            "genome_digest": GenomeDigest(D1),
            "genome_description": "human",
        },
        {
            "alias_name": GenomeAlias(D2),
            "genome_digest": GenomeDigest(D2),
            "genome_description": "no alias",
        },
    ]
    assert len(pulled) == 2
    puller.pull_multiple.assert_called_once()
