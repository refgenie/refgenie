"""
Unit tests for TransferManager's bulk-pull skeleton (``rgc.transfer``).

Every dependency is a MagicMock: no database, no server. These pin the shared
collect -> confirm -> pull_multiple path that pull_genomes and mirror use, and
how init_genomes and mirror hand each genome to ``init_from_remote``.
"""

import logging
from unittest.mock import MagicMock, call

import pytest

from refgenie.exceptions import PullFailedError
from refgenie.managers.transfer import TransferManager
from refgenie.managers.transfer.puller import BulkPullResult
from refgenie.models import GenomeAlias, GenomeDigest, first_usable_alias
from refgenie.plugins.events import EventSink

pytestmark = pytest.mark.unit

D1 = GenomeDigest("a" * 32)
D2 = GenomeDigest("b" * 32)


def _asset(group: str, digest: str = D1, size: int = 10) -> dict:
    return {"asset_group_name": group, "genome_digest": digest, "archive_size": size}


@pytest.fixture
def deps():
    servers = MagicMock()
    genome = MagicMock()
    puller = MagicMock()
    puller.pull_multiple.side_effect = lambda asset_list, **kwargs: BulkPullResult(
        pulled=list(asset_list)
    )
    return servers, genome, puller


@pytest.fixture
def transfer(deps):
    servers, genome, puller = deps
    return TransferManager(
        servers=servers,
        genome_manager=genome,
        puller=puller,
        events=EventSink(),
    )


def test_pull_genomes_filters_by_asset_group(transfer, deps):
    servers, _, puller = deps
    servers.list_assets_for_genome.return_value = [_asset("fasta"), _asset("bowtie2_index")]

    pulled = transfer.pull_genomes([D1], asset_group_name="fasta", force=True)

    assert pulled.pulled == [_asset("fasta")]
    assert puller.pull_multiple.call_args.kwargs["asset_list"] == [_asset("fasta")]


def test_pull_genomes_without_group_pulls_everything(transfer, deps):
    servers, _, puller = deps
    servers.list_assets_for_genome.return_value = [_asset("fasta"), _asset("bowtie2_index")]

    assert len(transfer.pull_genomes([D1], force=True).pulled) == 2
    puller.pull_multiple.assert_called_once()


def test_pull_genomes_all_genomes_expands_from_servers(transfer, deps):
    servers, genome, _ = deps
    servers.list_genomes.return_value = [
        {"aliases": ["hg38"], "genome_digest": D1},
        {"aliases": [], "genome_digest": D2},
    ]
    servers.list_assets_for_genome.side_effect = lambda d: [_asset("fasta", d)]

    pulled = transfer.pull_genomes(all_genomes=True, asset_group_name="fasta", force=True)

    # Every listed genome is used by digest; nothing is registered up front.
    genome.init_from_remote.assert_not_called()
    asked = [c.args[0] for c in servers.list_assets_for_genome.call_args_list]
    assert asked == [D1, D2]
    assert all(isinstance(d, GenomeDigest) for d in asked)
    assert len(pulled.pulled) == 2


def test_pull_genomes_registers_aliases_and_skips_unknown_ones(transfer, deps, caplog):
    servers, genome, _ = deps
    genome.init_from_remote.side_effect = [D1, None]
    servers.list_assets_for_genome.side_effect = lambda d: [_asset("fasta", d)]

    with caplog.at_level(logging.ERROR, logger="refgenie"):
        pulled = transfer.pull_genomes(
            [GenomeAlias("hg38"), GenomeAlias("nope"), D2], asset_group_name="fasta", force=True
        )

    assert genome.init_from_remote.call_args_list == [
        call(GenomeAlias("hg38")),
        call(GenomeAlias("nope")),
    ]
    assert [c.args[0] for c in servers.list_assets_for_genome.call_args_list] == [D1, D2]
    assert "Could not resolve genome: nope" in caplog.text
    assert len(pulled.pulled) == 2


def test_refused_confirmation_pulls_nothing(transfer, deps, caplog):
    servers, _, puller = deps
    servers.list_assets_for_genome.return_value = [_asset("fasta")]

    with caplog.at_level(logging.INFO, logger="refgenie"):
        assert transfer.pull_genomes([D1], confirm=lambda _msg: False) == BulkPullResult()

    puller.pull_multiple.assert_not_called()
    assert "Bulk pull cancelled" in caplog.text


def test_no_assets_warns_and_returns_empty(transfer, deps, caplog):
    servers, _, puller = deps
    servers.list_assets_for_genome.return_value = [_asset("bowtie2_index")]

    with caplog.at_level(logging.WARNING, logger="refgenie"):
        assert transfer.pull_genomes([D1], asset_group_name="fasta", force=True) == BulkPullResult()

    puller.pull_multiple.assert_not_called()
    assert "No remote assets named 'fasta' found for specified genomes" in caplog.text


def test_init_genomes_all_registers_aliasless_remote_genomes_by_digest(transfer, deps):
    """Every listed genome goes by digest; ``init_from_remote`` picks its alias."""
    servers, genome, _ = deps
    servers.list_genomes.return_value = [
        {"aliases": ["hg38"], "genome_digest": D1},
        {"aliases": [], "genome_digest": D2},
    ]
    genome.init_from_remote.side_effect = lambda g, **kwargs: g

    assert transfer.init_genomes(all_genomes=True) == [D1, D2]
    assert genome.init_from_remote.call_args_list == [
        call(D1, require_collection=True),
        call(D2, require_collection=True),
    ]


def test_init_genomes_tries_every_genome_then_raises_for_failures(transfer, deps):
    """A genome that cannot be registered does not stop the rest; failures raise (#299)."""
    _, genome, _ = deps
    genome.init_from_remote.side_effect = [D1, None, D2]

    with pytest.raises(PullFailedError, match="t7"):
        transfer.init_genomes([GenomeAlias("hg38"), GenomeAlias("t7"), D2])
    assert genome.init_from_remote.call_count == 3


def test_mirror_registers_each_genome_by_digest(transfer, deps):
    servers, genome, puller = deps
    servers.list_genomes.return_value = [
        {"aliases": ["hg38"], "genome_digest": D1, "description": "human"},
        {"aliases": [], "genome_digest": D2, "description": "no alias"},
    ]
    servers.list_assets_for_genome.side_effect = lambda d: [_asset("fasta", d)]

    pulled = transfer.mirror(force=True)

    assert genome.init_from_remote.call_args_list == [call(D1), call(D2)]
    assert all(isinstance(c.args[0], GenomeDigest) for c in genome.init_from_remote.call_args_list)
    assert len(pulled.pulled) == 2
    puller.pull_multiple.assert_called_once()


def test_first_usable_alias_skips_slash_and_empty_names():
    assert first_usable_alias(["a/b", "", "hg38", "mm10"]) == "hg38"
    assert first_usable_alias([]) is None
    assert first_usable_alias(["a/b"]) is None


@pytest.mark.parametrize("value", ["hg38", "a/b", "", "store::x", None, 3])
def test_genome_alias_is_valid_agrees_with_the_constructor(value):
    if GenomeAlias.is_valid(value):
        GenomeAlias(value)
    else:
        with pytest.raises(ValueError):
            GenomeAlias(value)
