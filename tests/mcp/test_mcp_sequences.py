"""MCP sequence tools over real stores."""

from pathlib import Path

import pytest

from refgenie import Refgenie
from refgenie.db.tables import StoreType
from tests.support.factories import build_on_disk_store, make_engine

try:
    from mcp.server.mcpserver.exceptions import ToolError

    from refgenie.mcp import tools
except ImportError:
    tools = None

pytestmark = [
    pytest.mark.component,
    pytest.mark.skipif(tools is None, reason="requires the 'mcp' extras"),
]

FASTA = Path(__file__).resolve().parents[1] / "data" / "refget_test_fasta"


@pytest.fixture
def two_genomes(tmp_path):
    """A local Refgenie holding the base and different_names genomes."""
    rg = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    rg.database.init(genome_folder=tmp_path / "genomes")
    for name in ("base", "different_names"):
        rg.genome.initialize_genome(
            fasta_file_path=FASTA / f"{name}.fa", alias_names=[name], description=name
        )
    tools.set_refgenie(rg)
    yield rg
    tools.set_refgenie(None)


def test_compare_genomes_summarizes_and_maps_names(two_genomes):
    result = tools.compare_genomes(alias_a="base", alias_b="different_names")
    assert result.name_map == {"chrX": "X", "chr1": "1", "chr2": "2"}
    assert result.summary.startswith("Same 3 sequences in the same order; all 3 names differ")
    assert "array_elements" in result.comparison


def test_get_sequence_reads_a_region(two_genomes):
    assert tools.get_sequence(alias="base", name="chrX", start=0, end=4) == "TTGG"
    with pytest.raises(ToolError, match="nope"):
        tools.get_sequence(alias="base", name="nope", start=0, end=4)


def test_sequence_tools_work_in_server_mode(tmp_path):
    from refgenie.cli.commands.store import _sync_one_store

    path = tmp_path / "store"
    digest = build_on_disk_store(path, FASTA / "base.fa", "base")
    rg = Refgenie(database_engine=make_engine(), suppress_migrations=True, server_mode=True)
    rg.database.create_tables()
    rg.store.apply_declaration(
        [
            {
                "name": "jungle",
                "url": str(path),
                "type": StoreType.on_disk.value,
                "priority": 10,
                "enabled": True,
            }
        ]
    )
    _sync_one_store(rg, rg.store.get("jungle"), {"jungle": 10}, 100)
    tools.set_refgenie(rg, public_url="https://rg.example")
    try:
        assert tools.get_sequence(digest=digest, name="chrX", start=0, end=4) == "TTGG"
        page = tools.get_genome_sequences(digest=digest)
        assert page.n_sequences == 3 and page.results[0].name == "chrX"
    finally:
        tools.set_refgenie(None)
