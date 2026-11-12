"""The MCP tools (``refgenie.mcp.tools``) and the import guard behind them.

The guard tests need the ``mcp`` extra absent, so the import below is a plain
try/except rather than a module-level ``importorskip``.
"""

import asyncio
import importlib
import sys
from pathlib import Path

import pytest
from sqlmodel import Session

from refgenie.db.tables import Asset, AssetGroup, AssetName
from tests.support.factories import register_fasta, seed_search_genomes, seed_species

try:
    from mcp.server.mcpserver.exceptions import ToolError

    from refgenie.mcp import tools
except ImportError:  # the 'mcp' extras are not installed
    tools = None

requires_mcp = pytest.mark.skipif(tools is None, reason="requires the 'mcp' extras")

A = "a" * 32
B = "b" * 32


# ---------------------------------------------------------------------------
# The tool implementations (need the extras)
# ---------------------------------------------------------------------------


@pytest.fixture
def mcp_refgenie(refgenie_minimal):
    """A Refgenie with two genomes whose searchable fields are all distinct."""
    seed_search_genomes(refgenie_minimal)

    tools.set_refgenie(refgenie_minimal)
    yield refgenie_minimal
    tools.set_refgenie(None)


@pytest.fixture
def chain(mcp_refgenie):
    """Three fasta assets of genome hg38 linked ``top -> mid -> leaf``.

    Rows only, no files: enough for the catalog reads. Returns the digests.
    """
    register_fasta(mcp_refgenie)
    asset_class = mcp_refgenie.asset_class.get("fasta")
    digests = {role: role[0] * 64 for role in ("top", "mid", "leaf")}
    with Session(mcp_refgenie.database_engine) as session:
        group = AssetGroup(name="fasta", genome_digest=A, asset_class_id=asset_class.id)
        session.add(group)
        session.flush()
        for role, digest in digests.items():
            session.add(
                Asset(
                    digest=digest,
                    name=role,
                    asset_group_id=group.id,
                    path=f"data/{A}/fasta/{digest}",
                    size=4,
                )
            )
            session.add(AssetName(name=role, asset_group_id=group.id, asset_digest=digest))
        session.flush()
        top, mid, leaf = (session.get(Asset, digests[r]) for r in ("top", "mid", "leaf"))
        mid.parents.append(top)
        leaf.parents.append(mid)
        session.commit()
    return digests


def _search(query, **kwargs):
    return tools.search_genomes(query, **kwargs).results


@requires_mcp
@pytest.mark.parametrize("kwargs", [{"alias": "hg38"}, {"digest": A}], ids=["alias", "digest"])
def test_get_genome_by_alias_or_by_digest(mcp_refgenie, kwargs):
    """get_genome finds the same genome by its alias parameter or its digest parameter."""
    result = tools.get_genome(**kwargs)
    assert result.digest == A
    assert result.aliases == ["hg38"]
    assert result.asset_count == 0


@requires_mcp
@pytest.mark.parametrize("kwargs", [{}, {"alias": "hg38", "digest": A}], ids=["neither", "both"])
def test_get_genome_needs_exactly_one_identifier(mcp_refgenie, kwargs):
    """Neither or both of alias and digest is an error that says which to give."""
    with pytest.raises(ToolError, match="exactly one of alias"):
        tools.get_genome(**kwargs)


@requires_mcp
def test_a_digest_is_not_an_alias(mcp_refgenie):
    """A digest passed as an alias is looked up as an alias, is not one, and says what to try."""
    with pytest.raises(ToolError, match="not found.*search_genomes"):
        tools.get_genome(alias=A)


@requires_mcp
def test_get_genome_unknown_digest_raises(mcp_refgenie):
    """A well-formed digest of no genome here names the genome and the alias= escape."""
    with pytest.raises(ToolError, match="not found.*pass alias="):
        tools.get_genome(digest="c" * 32)


@requires_mcp
def test_a_name_passed_as_a_digest_is_a_readable_error(mcp_refgenie):
    with pytest.raises(ToolError, match="'hg38' is not a genome digest.*pass alias="):
        tools.get_genome(digest="hg38")


@requires_mcp
def test_list_assets_and_compare_take_explicit_identifiers(mcp_refgenie):
    """list_assets filters by alias or digest; compare_genomes names each side."""
    assert tools.list_assets(alias="hg38").results == []
    assert tools.list_assets(digest=B).results == []
    with pytest.raises(ToolError, match="exactly one of alias_b"):
        tools.compare_genomes(alias_a="hg38")


@requires_mcp
@pytest.mark.parametrize(
    "query",
    ["sapiens", "Peculiar", "UCSC", "GCA_000001405", "hg38", "aaaa"],
    ids=[
        "species",
        "description",
        "assembly_source",
        "assembly_accession",
        "alias",
        "digest",
    ],
)
def test_search_matches_every_searchable_field(mcp_refgenie, query):
    """Each searchable field finds the human genome and only the human genome."""
    assert [g.digest for g in _search(query)] == [A]


@requires_mcp
@pytest.mark.parametrize("query", ["SAPIENS", "sapiens", "SaPiEnS"])
def test_search_is_case_insensitive(mcp_refgenie, query):
    """Matching ignores case."""
    assert [g.digest for g in _search(query)] == [A]


@requires_mcp
def test_search_returns_slim_genome_rows(mcp_refgenie):
    """A row is the small shape list_genomes uses, not the full record."""
    (result,) = _search("hg38")
    assert result.model_dump() == {
        "digest": A,
        "aliases": ["hg38"],
        "species_name": "Homo sapiens",
        "common_name": None,
        "asset_count": 0,
    }


@requires_mcp
def test_search_can_match_multiple_genomes(mcp_refgenie):
    """A term shared by both genomes returns both, once each."""
    assert sorted(g.digest for g in _search("m")) == [A, B]


@requires_mcp
def test_search_with_no_match_returns_empty_with_a_hint(mcp_refgenie):
    page = tools.search_genomes("no-such-genome")
    assert page.results == []
    assert "names, not meaning" in page.hint
    assert tools.search_genomes("hg38").hint is None


@requires_mcp
def test_search_finds_genome_synonyms(mcp_refgenie):
    """A synonym from FHR genomeSynonym is searched, and get_genome returns it."""
    mcp_refgenie.genome.apply_fhr_columns(B, {"genomeSynonym": ["GRCm39"]})
    assert [g.digest for g in _search("GRCm")] == [B]
    assert tools.get_genome(alias="mm39").synonyms == ["GRCm39"]
    assert tools.get_genome(alias="hg38").synonyms is None


@requires_mcp
def test_search_finds_species_list_names(mcp_refgenie):
    """A name from the species list finds the genome; get_genome returns the species' names."""
    seed_species(mcp_refgenie, 10090, "Mus musculus", "house mouse", ["field mouse"])
    assert [g.digest for g in _search("field mouse")] == [B]
    assert tools.get_genome(alias="mm39").species.synonyms == ["field mouse"]
    assert tools.get_genome(alias="hg38").species is None


@requires_mcp
def test_list_rows_show_the_species_common_name(mcp_refgenie):
    """List rows show the common name of the genome's species row, if it has one."""
    seed_species(mcp_refgenie, 10090, "Mus musculus", "house mouse")
    rows = {g.digest: g for g in tools.list_genomes().results}
    assert rows[B].common_name == "house mouse"
    assert rows[A].common_name is None


@requires_mcp
def test_list_genomes_filters_by_taxon_and_assets(mcp_refgenie, chain):
    assert [g.digest for g in tools.list_genomes(taxon_id=9606).results] == [A]
    assert [g.digest for g in tools.list_genomes(has_assets=True).results] == [A]
    assert [g.digest for g in tools.list_genomes(has_assets=False).results] == [B]
    # Genomes with assets come first.
    assert [g.digest for g in tools.list_genomes().results] == [A, B]


@requires_mcp
def test_token_round_trip_walks_every_row_once(mcp_refgenie):
    first = tools.list_genomes(page_size=1)
    assert first.pagination.page_size == 1
    assert first.pagination.total == 2
    second = tools.list_genomes(page_size=1, token=first.pagination.next_page_token)
    assert second.pagination.next_page_token is None
    assert {first.results[0].digest, second.results[0].digest} == {A, B}


@requires_mcp
def test_token_from_other_filters_is_refused(mcp_refgenie):
    token = tools.list_genomes(page_size=1).pagination.next_page_token
    with pytest.raises(ToolError, match="does not match this query"):
        tools.list_genomes(page_size=1, taxon_id=9606, token=token)
    with pytest.raises(ToolError, match="does not match this query"):
        tools.search_genomes("m", page_size=1, token=token)


@requires_mcp
@pytest.mark.parametrize("token", ["garbage", "e30", "!!!"])
def test_garbage_token_is_refused(mcp_refgenie, token):
    with pytest.raises(ToolError, match="not one this server issued"):
        tools.list_genomes(token=token)


@requires_mcp
def test_page_size_is_clamped_and_echoed(mcp_refgenie):
    assert tools.list_genomes(page_size=500).pagination.page_size == 100
    assert tools.list_genomes(page_size=0).pagination.page_size == 20
    assert tools.list_genomes().pagination.page_size == 20


@requires_mcp
def test_get_asset_loads_parents_and_children(mcp_refgenie, chain):
    """The asset in the middle of a chain reports both relatives."""
    asset = tools.get_asset(chain["mid"])
    assert [p.digest for p in asset.parents] == [chain["top"]]
    assert [c.digest for c in asset.children] == [chain["leaf"]]
    assert asset.asset_class == "fasta"
    assert asset.names == ["mid"]
    assert asset.pull_command == "refgenie pull hg38/fasta:mid"
    # stdio: no archive URL, and a chain asset with no seek keys has no files.
    assert asset.archive_url is None
    assert asset.files == {}
    assert asset.description == mcp_refgenie.asset_class.get("fasta").description


@requires_mcp
def test_get_asset_on_a_server_gives_download_urls(mcp_refgenie, chain):
    tools.set_refgenie(mcp_refgenie, public_url="https://rg.example")
    asset = tools.get_asset(chain["mid"])
    assert asset.archive_url == f"https://rg.example/v4/archives/{chain['mid']}/download"


def _add_seek_keys(rgc, digest, **keys):
    """Give an asset seek keys: ``name=(value, type)``."""
    from refgenie.db.tables import SeekKey, SeekKeyType

    with Session(rgc.database_engine) as session:
        asset = session.get(Asset, digest)
        folder = Path(rgc.genome_folder) / asset.path
        folder.mkdir(parents=True, exist_ok=True)
        for name, (value, kind) in keys.items():
            (folder / value).write_text("x")
            session.add(
                SeekKey(
                    name=name,
                    value=value,
                    description=None,
                    type=SeekKeyType(kind),
                    asset=asset,
                )
            )
        session.commit()


@requires_mcp
def test_pull_command_uses_the_shortest_plain_alias(mcp_refgenie, chain, monkeypatch):
    mcp_refgenie.alias.add("GRCh38-igenomes-ncbi", A)
    mcp_refgenie.alias.add("hs38", A)
    assert tools.get_asset(chain["mid"]).pull_command == "refgenie pull hg38/fasta:mid"
    # An alias with "/" is no use in a command; with no plain alias, use the digest.
    monkeypatch.setattr(mcp_refgenie.alias, "get_for_genomes", lambda digests: {A: ["a/b"]})
    expected = f"refgenie pull fasta:mid --genome-digest {A}"
    assert tools.get_asset(chain["mid"]).pull_command == expected


@requires_mcp
def test_files_give_seek_paths_locally_and_urls_on_a_server(mcp_refgenie, chain):
    mid = chain["mid"]
    _add_seek_keys(mcp_refgenie, mid, fa=("x y.fa", "file"), idx=("x", "prefix"))
    files = tools.get_asset(mid).files
    for key in ("fa", "idx"):
        assert files[key] == mcp_refgenie.asset.seek(A, "fasta", "mid", key)
    tools.set_refgenie(mcp_refgenie, public_url="https://rg.example")
    assert tools.get_asset(mid).files == {
        "fa": f"https://rg.example/v4/assets/{mid}/files/x%20y.fa"
    }


@requires_mcp
def test_hints_for_empty_results(mcp_refgenie):
    hint = tools.list_assets(alias="mm39").hint
    assert f"template=asset_request.yml&genome=mm39&genome_digest={B}" in hint
    assert "genome_proposal.yml&genome=zzz" in tools.search_genomes("zzz").hint
    assert "4932" in tools.list_genomes(taxon_id=559292).hint


@requires_mcp
def test_catalog_summary_and_instructions_name_the_server(mcp_refgenie):
    url = mcp_refgenie.database_engine.url.render_as_string(hide_password=True)
    assert url in tools.catalog_summary().server
    assert url in tools.mcp.instructions
    tools.set_refgenie(mcp_refgenie, public_url="https://rg.example")
    assert "https://rg.example" in tools.catalog_summary().server
    assert tools.mcp.instructions.startswith("refgenie server at https://rg.example.")
    tools.set_refgenie(None)
    assert tools.mcp.instructions == tools.INSTRUCTIONS


@requires_mcp
def test_get_sequence_refuses_more_than_the_limit(mcp_refgenie):
    with pytest.raises(ToolError, match="At most 100,000 bases"):
        tools.get_sequence("chr1", 0, 100_001, alias="hg38")


@requires_mcp
def test_get_genome_sequences_pages_and_filters(mcp_refgenie, monkeypatch):
    level2 = {
        "names": ["chr1", "chr2", "chrM"],
        "lengths": [10, 20, 5],
        "sequences": ["SQ.1", "SQ.2", "SQ.m"],
    }
    monkeypatch.setattr(mcp_refgenie.genome, "level2", lambda digest: level2)
    page = tools.get_genome_sequences(alias="hg38", name_contains="CHR")
    assert [r.name for r in page.results] == ["chr1", "chr2", "chrM"]
    assert (page.n_sequences, page.total_length) == (3, 35)
    assert [
        r.name for r in tools.get_genome_sequences(alias="hg38", name_contains="m").results
    ] == ["chrM"]
    first = tools.get_genome_sequences(alias="hg38", page_size=1)
    assert first.results[0].digest == "SQ.1" and first.pagination.total == 3
    second = tools.get_genome_sequences(
        alias="hg38", page_size=1, token=first.pagination.next_page_token
    )
    assert second.results[0].name == "chr2"
    with pytest.raises(ToolError, match="does not match"):
        tools.get_genome_sequences(
            alias="hg38", name_contains="x", token=first.pagination.next_page_token
        )


@requires_mcp
def test_get_asset_unknown_digest_hints_at_list_assets(mcp_refgenie):
    with pytest.raises(ToolError, match="not found.*list_assets"):
        tools.get_asset("z" * 64)


@requires_mcp
def test_lookup_digest_says_what_a_digest_is(mcp_refgenie, chain):
    assert tools.lookup_digest(A).type == "genome"
    found = tools.lookup_digest(chain["leaf"])
    assert (found.type, found.asset.name) == ("asset", "leaf")
    assert tools.lookup_digest("nope").type == "not_found"


@requires_mcp
def test_recipes_and_asset_classes_are_slim_and_paged(mcp_refgenie):
    register_fasta(mcp_refgenie)
    (recipe,) = tools.list_recipes().results
    assert (recipe.name, recipe.output_asset_class) == ("fasta", "fasta")
    assert tools.list_recipes(output_asset_class="nothing").results == []
    (asset_class,) = tools.list_asset_classes().results
    assert asset_class.name == "fasta"
    assert asset_class.seek_keys


@requires_mcp
def test_get_recipe_defaults_to_the_latest_version(mcp_refgenie):
    register_fasta(mcp_refgenie)
    recipe = tools.get_recipe("fasta")
    assert recipe.version == "0.1.0"
    assert recipe.command_templates
    assert recipe.next == "refgenie build <genome>/fasta"
    with pytest.raises(ToolError, match="list_recipes"):
        tools.get_recipe("no-such-recipe")


@requires_mcp
def test_catalog_summary_counts_species_with_assets(mcp_refgenie, chain):
    summary = tools.catalog_summary()
    assert (summary.genomes, summary.asset_groups, summary.assets) == (2, 1, 3)
    assert [(s.species, s.assets) for s in summary.species_with_assets] == [("Homo sapiens", 3)]


@requires_mcp
def test_errors_reach_the_agent_through_the_protocol(mcp_refgenie):
    """Through the MCP server, a refgenie error keeps its message and hint.

    Any exception other than ``ToolError`` would reach the agent as a bare
    "Error executing tool get_genome", with the message withheld.
    """
    with pytest.raises(ToolError) as excinfo:
        asyncio.run(tools.mcp.call_tool("get_genome", {"alias": "hg39"}))
    assert "Genome alias 'hg39' not found" in str(excinfo.value)
    assert "search_genomes" in str(excinfo.value)


def _call(name, arguments):
    return asyncio.run(tools.mcp.call_tool(name, arguments))


@requires_mcp
def test_an_unknown_argument_is_refused_not_dropped(mcp_refgenie):
    """mcp would drop genome= and return every asset; the agent would believe
    it got hg38's."""
    with pytest.raises(ToolError) as excinfo:
        _call("list_assets", {"genome": "hg38"})
    message = str(excinfo.value)
    assert "list_assets has no argument 'genome'" in message
    assert "alias, digest, asset_class, page_size, token" in message
    assert "alias= (like hg38)" in message


@requires_mcp
def test_a_misspelled_argument_gets_a_suggestion(mcp_refgenie):
    with pytest.raises(ToolError, match="Did you mean 'page_size'\\?"):
        _call("list_genomes", {"page_sise": 5})


@requires_mcp
def test_known_arguments_still_work_through_the_protocol(mcp_refgenie):
    assert not _call("list_genomes", {"page_size": 1, "taxon_id": 9606}).is_error


@requires_mcp
def test_every_tool_has_its_parameters_recorded(mcp_refgenie):
    for tool in asyncio.run(tools.mcp.list_tools()):
        assert tool.name in tools._PARAMS
        assert set(tool.input_schema.get("properties", {})) == set(tools._PARAMS[tool.name])


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
        stdio.main([])
    assert excinfo.value.code == 1
    _assert_friendly(capsys.readouterr().err)


def test_refgenie_mcp_help_works_without_the_extra(without_mcp, capsys):
    stdio = importlib.import_module("refgenie.mcp.stdio")
    with pytest.raises(SystemExit) as excinfo:
        stdio.main(["--help"])
    assert excinfo.value.code == 0
    assert "REFGENIE_HOME_PATH" in capsys.readouterr().out


@requires_mcp
@pytest.mark.parametrize(
    "doc", ["docs/site/SKILL.md", "frontend/public/SKILL.md", "docs/refgenie/mcp.md"]
)
def test_doc_tool_tables_match_registered_tools(doc):
    import re

    root = Path(__file__).resolve().parents[2]
    lines = (root / doc).read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("| Tool |"))
    documented = set()
    for line in lines[start + 2 :]:
        if not line.startswith("|"):
            break
        documented.update(re.findall(r"`(\w+)`", line.split("|")[1]))
    assert documented == {t.name for t in asyncio.run(tools.mcp.list_tools())}
