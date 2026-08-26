"""
Tests for refgenie's MCP surface: the tool implementations
(`refgenie.mcp.tools`) and the import guard behind them.

`search_genomes` runs its substring match in SQL for scalar columns, plus a
mode-selected alias-manager lookup for aliases; it must cover every
searchable field, case-insensitively. Aliases are added through
`rgc.alias.add()` (the mode-agnostic manager API) rather than inserted
directly into the SQL `Alias` table, since `refgenie_minimal` runs in local
(store-backed) mode where the SQL alias table is never read.

The `refgenie-mcp` console script is installed by every base `pip install
refgenie`, but the `mcp` package it needs only ships in the `[mcp]` extra.
Importing the entry point without that extra must produce the same friendly,
actionable ImportError that `refgenie serve` and `refgenie dash` produce, not a
bare ModuleNotFoundError. `mcp` is generally installed in the dev environment,
so absence is simulated with a meta path finder that blocks it.

The tool tests need the extras; the guard tests need them *absent*. So the
import below is a plain try/except rather than a module-level `importorskip`,
which would skip the guard tests -- the ones that matter without the extras --
along with everything else.
"""

import importlib
import json
import sys

import pytest
from sqlmodel import Session

from refgenie.db.tables import Genome

try:
    from refgenie.mcp import tools
except ImportError:  # the 'mcp' extras are not installed
    tools = None

requires_mcp = pytest.mark.skipif(tools is None, reason="requires the 'mcp' extras")


# ---------------------------------------------------------------------------
# The tool implementations (need the extras)
# ---------------------------------------------------------------------------


@pytest.fixture
def mcp_refgenie(refgenie_minimal):
    """A Refgenie with two genomes whose searchable fields are all distinct."""
    with Session(refgenie_minimal.database_engine) as session:
        session.add(
            Genome(
                digest="a" * 32,
                description="A Peculiar Description",
                species_name="Homo sapiens",
                common_name="Human",
                taxon_id=9606,
                assembly_source="UCSC",
                assembly_accession="GCA_000001405.15",
            )
        )
        session.add(
            Genome(
                digest="b" * 32,
                description="unrelated",
                species_name="Mus musculus",
                common_name="Mouse",
                taxon_id=10090,
                assembly_source="Ensembl",
                assembly_accession="GCA_000001635.9",
            )
        )
        session.commit()

    # Aliases go through the mode-selected manager, not the SQL Alias table:
    # refgenie_minimal is local (store-backed) mode.
    refgenie_minimal.alias.add("hg38", "a" * 32)
    refgenie_minimal.alias.add("mm39", "b" * 32)

    tools.set_refgenie(refgenie_minimal)
    yield refgenie_minimal
    tools.set_refgenie(None)


def _search(query):
    return json.loads(tools.search_genomes(query))


@requires_mcp
@pytest.mark.parametrize(
    "query",
    ["sapiens", "Human", "Peculiar", "UCSC", "GCA_000001405", "hg38"],
    ids=[
        "species_name",
        "common_name",
        "description",
        "assembly_source",
        "assembly_accession",
        "alias",
    ],
)
def test_search_matches_every_searchable_field(mcp_refgenie, query):
    """Each searchable field finds the human genome and only the human genome."""
    results = _search(query)
    assert [g["digest"] for g in results] == ["a" * 32]


@requires_mcp
@pytest.mark.parametrize("query", ["SAPIENS", "sapiens", "SaPiEnS"])
def test_search_is_case_insensitive(mcp_refgenie, query):
    """Matching ignores case."""
    assert [g["digest"] for g in _search(query)] == ["a" * 32]


@requires_mcp
def test_search_returns_full_genome_payload(mcp_refgenie):
    """The result shape matches list_genomes' entries exactly."""
    (result,) = _search("hg38")
    assert result == {
        "digest": "a" * 32,
        "aliases": ["hg38"],
        "species_name": "Homo sapiens",
        "common_name": "Human",
        "taxon_id": 9606,
        "assembly_source": "UCSC",
        "assembly_accession": "GCA_000001405.15",
        "description": "A Peculiar Description",
    }


@requires_mcp
def test_search_can_match_multiple_genomes(mcp_refgenie):
    """A term shared by both genomes returns both, once each."""
    assert sorted(g["digest"] for g in _search("m")) == ["a" * 32, "b" * 32]


@requires_mcp
def test_search_with_no_match_returns_empty(mcp_refgenie):
    assert _search("no-such-genome") == []


@requires_mcp
def test_get_genome_with_no_asset_groups(mcp_refgenie):
    """A genome with zero asset groups must not raise DetachedInstanceError."""
    result = json.loads(tools.get_genome("hg38"))
    assert result["digest"] == "a" * 32
    assert result["asset_groups"] == []
    assert result["aliases"] == ["hg38"]


# ---------------------------------------------------------------------------
# The import guard (needs the extras to be ABSENT -- simulated)
# ---------------------------------------------------------------------------


class _BlockModule:
    """Meta path finder that makes a top-level module (and submodules) unimportable."""

    def __init__(self, name):
        self.name = name

    def find_spec(self, fullname, path=None, target=None):
        if fullname == self.name or fullname.startswith(self.name + "."):
            raise ModuleNotFoundError(f"No module named '{fullname}'", name=fullname)
        return None


@pytest.fixture
def without_mcp():
    """Make `import mcp` fail, with sys.modules restored afterwards."""
    purged = {
        name: mod
        for name, mod in list(sys.modules.items())
        if name == "mcp" or name.startswith("mcp.") or name.startswith("refgenie.mcp")
    }
    for name in purged:
        del sys.modules[name]

    finder = _BlockModule("mcp")
    sys.meta_path.insert(0, finder)
    try:
        yield
    finally:
        sys.meta_path.remove(finder)
        for name in [
            n
            for n in list(sys.modules)
            if n == "mcp" or n.startswith("mcp.") or n.startswith("refgenie.mcp")
        ]:
            del sys.modules[name]
        sys.modules.update(purged)


def _assert_friendly(message):
    assert "mcp" in message
    assert "extras" in message
    assert "install" in message.lower()
    # A bare ModuleNotFoundError leaking through is the bug we're guarding against.
    assert "No module named" not in message


def test_mcp_tools_import_gives_friendly_error(without_mcp):
    """refgenie.mcp.tools must explain which extra is missing."""
    with pytest.raises(ImportError) as excinfo:
        importlib.import_module("refgenie.mcp.tools")
    _assert_friendly(str(excinfo.value))


def test_mcp_stdio_entry_point_exits_with_friendly_error(without_mcp, capsys):
    """`refgenie-mcp` must exit non-zero with an actionable message, not a traceback."""
    stdio = importlib.import_module("refgenie.mcp.stdio")
    with pytest.raises(SystemExit) as excinfo:
        stdio.main()
    assert excinfo.value.code == 1
    _assert_friendly(capsys.readouterr().err)
