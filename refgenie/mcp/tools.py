"""
MCP tool definitions for Refgenie.

All tools are read-only queries against the Refgenie object, through the same
manager reads the HTTP API uses; they present the results as the slim models
in ``schemas.py``. Shared between stdio (local) and Streamable HTTP (server)
transports.
"""

import difflib
import functools
import inspect
from urllib.parse import quote

try:
    from mcp.server.mcpserver import MCPServer
    from mcp.server.mcpserver.exceptions import ToolError
except ImportError as e:
    raise ImportError(
        "The 'mcp' extras are not installed. Please install refgenie with the 'mcp' extras to use the MCP server."
    ) from e

from refgenie import __version__ as _refgenie_version
from refgenie.db.tables import Asset, SeekKeyType, is_path_type
from refgenie.exceptions import (
    MissingAliasError,
    MissingAssetError,
    MissingGenomeError,
    MissingRecipeError,
    RefgenieError,
)
from refgenie.mcp import paging
from refgenie.mcp.paging import ResultPage
from refgenie.mcp.schemas import (
    AssetClassRow,
    AssetDetail,
    AssetRow,
    CatalogSummary,
    DigestLookup,
    GenomeComparison,
    GenomeDetail,
    GenomeRow,
    RecipeDetail,
    RecipeRow,
    SequencePage,
    SequenceRow,
    SpeciesCount,
    SpeciesNames,
)
from refgenie.models import GenomeAlias, GenomeDigest, Page, display_alias
from refgenie.utils.search import Search
from refgenie.utils.seqcol import compare_level2, summarize_comparison

INSTRUCTIONS = """\
refgenie is a catalog of reference genomes and the assets built from them (FASTA files, aligner indexes, annotations).
A genome has a 32-character digest and any number of aliases (names like hg38). Tools take either alias= or digest=, never a guess between them.
Most genomes have no assets; pass has_assets=true to list only genomes with files.
Which tool for which job:
- Find a genome: search_genomes or list_genomes (has_assets=true for genomes with files).
- Same genome? compare_genomes; summary says what is shared, name_map maps renamed sequences.
- Chromosome names and sizes: get_genome_sequences.
- Bases of a region: get_sequence, 0-based half-open like BED, at most 100,000 bases.
- Files: list_assets(alias=...), then get_asset(digest) for files, the archive URL and the pull command.
search_genomes matches every name of a species, so other common names work: brown rat finds Rattus norvegicus even where its common name is listed as Norway rat.
taxon_id= filters genomes by NCBI taxon id (9606 is human).
List tools return 20 rows per page (page_size up to 100); pass pagination.next_page_token back as token for the next page.
Every tool is read-only.
"""

#: Each tool's parameter names, to refuse arguments a tool does not take.
_PARAMS: dict[str, tuple[str, ...]] = {}


def _unknown_arguments(tool: str, unknown: list[str], params: tuple[str, ...]) -> str:
    """The error for arguments ``tool`` does not take: what was wrong, what is
    allowed, and what to do instead."""
    if not params:
        return f"{tool} takes no arguments."
    message = [
        f"{tool} has no argument {', '.join(repr(a) for a in unknown)}. "
        f"Its arguments: {', '.join(params)}."
    ]
    for arg in unknown:
        if close := difflib.get_close_matches(arg, params, n=1):
            message.append(f"Did you mean '{close[0]}'?")
    if "alias" in params:
        message.append("Name a genome with alias= (like hg38) or digest= (32 characters).")
    return " ".join(message)


class _RefgenieMCP(MCPServer):
    """An MCPServer that refuses arguments a tool does not take.

    mcp drops unknown arguments without a word, so list_assets(genome="hg38")
    would return every asset.
    """

    async def call_tool(self, name, arguments, context=None):
        params = _PARAMS.get(name)
        if params is not None and (unknown := sorted(set(arguments) - set(params))):
            raise ToolError(_unknown_arguments(name, unknown, params))
        return await super().call_tool(name, arguments, context)


mcp = _RefgenieMCP("refgenie", version=_refgenie_version, instructions=INSTRUCTIONS)

#: Longest region get_sequence returns; BRCA1, about 81 kb, fits.
MAX_SEQUENCE_REGION = 100_000

_REGISTRY_ISSUES = "https://github.com/refgenie/refgenie-registry/issues/new"

_rgc = None
#: None in local mode; the server's base URL in server mode.
_public_url: str | None = None

#: What to try next, per error type, appended to the error the agent reads.
_HINTS: tuple[tuple[type[Exception], str], ...] = (
    (MissingAliasError, "Find genomes with search_genomes(query) or list_genomes()."),
    (MissingGenomeError, "A digest is 32 characters; for a name like hg38 pass alias=."),
    (MissingAssetError, "list_assets(alias=...) shows a genome's assets."),
    (MissingRecipeError, "list_recipes() shows recipe names."),
)

_NO_MATCH_HINT = (
    "No match. Search matches names, not meaning: try the scientific name (Rattus norvegicus) "
    "or another common name (brown rat), or taxon_id with the NCBI taxon id (10116). "
    "has_assets=true shows only genomes with files. If refgenie does not have the genome, "
    f"request it: {_REGISTRY_ISSUES}?template=genome_proposal.yml&genome={{query}}."
)

_TAXON_HINT = (
    "No genome is tagged with NCBI taxon {taxon_id}. Genomes carry their species' taxon id: "
    "for a strain, use the species id (S. cerevisiae S288C is 4932, not 559292), "
    "or search_genomes with the species name."
)

_NO_ASSETS_HINT = (
    "No assets for {label} here. On the public server (https://api.refgenie.org/mcp) "
    "list_assets shows what you can pull. To ask for builds, open "
    f"{_REGISTRY_ISSUES}?template=asset_request.yml&genome={{alias}}&genome_digest={{digest}}."
)


def _get_refgenie():
    """The Refgenie the tools read, set by ``set_refgenie``."""
    return _rgc


def _server_line() -> str:
    """One sentence saying which server this is."""
    if _public_url is not None:
        return f"refgenie server at {_public_url}."
    url = _rgc.database_engine.url.render_as_string(hide_password=True)
    return (
        f"Local refgenie, database {url}. It lists only what is on this machine; "
        "refgenie's full catalog is at https://api.refgenie.org/mcp."
    )


def set_refgenie(rgc, *, public_url: str | None = None):
    """Set the Refgenie the tools read, and the mode.

    Args:
        rgc: The Refgenie instance the tools read.
        public_url: None for local mode. For a server, its base URL, which
            ``get_asset`` uses for download URLs.
    """
    global _rgc, _public_url
    _rgc, _public_url = rgc, public_url
    # The one private attribute used; test_initialize_names_the_server guards it.
    mcp._lowlevel_server.instructions = (
        f"{_server_line()}\n{INSTRUCTIONS}" if rgc is not None else INSTRUCTIONS
    )


def _tool(fn):
    """Register ``fn`` as a tool whose expected failures reach the agent.

    mcp reports any exception other than ``ToolError`` as a bare "Error
    executing tool <name>", so a refgenie error or a bad argument is re-raised
    as a ``ToolError`` carrying its message and a hint at the next call.
    """
    _PARAMS[fn.__name__] = tuple(inspect.signature(fn).parameters)

    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except (RefgenieError, ValueError) as e:
            hint = next((h for t, h in _HINTS if isinstance(e, t)), "")
            raise ToolError(f"{e} {hint}".strip()) from e

    return mcp.tool()(wrapper)


def _genome_digest(alias: str | None, digest: str | None, label: str = "") -> GenomeDigest:
    """The digest of the genome a tool call names, by alias or by digest.

    Exactly one of ``alias`` and ``digest`` must be given. Nothing is guessed:
    an alias is resolved through the alias manager, a digest is used as is.
    """
    if (alias is None) == (digest is None):
        raise ValueError(
            f"Give exactly one of alias{label} (a genome name like 'hg38') "
            f"or digest{label} (a 32-character genome digest)."
        )
    if digest is not None:
        try:
            return GenomeDigest(digest)
        except ValueError as e:
            raise ValueError(f"{e}. For a name like hg38 pass alias{label}=.") from None
    return _get_refgenie().alias.resolve(GenomeAlias(alias))


def _asset_detail(asset: Asset) -> AssetDetail:
    """An asset's full record, with how to get its files from here."""
    rgc = _get_refgenie()
    aliases = rgc.alias.get_for_genomes([asset.genome_digest]).get(asset.genome_digest, [])
    group, genome = asset.asset_group_name, asset.genome_digest
    alias = display_alias(aliases)
    target = (
        f"{alias}/{group}:{asset.name}"
        if alias
        else f"{group}:{asset.name} --genome-digest {genome}"
    )
    # Over stdio the asset is already here, so say how to find it, not how to pull it.
    local = _public_url is None
    archive_url = None
    files: dict[str, str] = {}
    if local:
        files = {
            sk.name: rgc.asset.seek(genome, group, asset.name, sk.name)
            for sk in asset.seek_keys
            if is_path_type(sk.type)
        }
    else:
        if "archive" in asset.serving_modes:
            archive_url = f"{_public_url}/v4/archives/{asset.digest}/download"
        if "file" in asset.serving_modes:
            files = {
                sk.name: f"{_public_url}/v4/assets/{asset.digest}/files/{quote(sk.value)}"
                for sk in asset.seek_keys
                if sk.type == SeekKeyType.file
            }
    asset_class = asset.asset_group.asset_class if asset.asset_group else None
    return AssetDetail(
        **AssetRow.of(asset).model_dump(),
        asset_group=group,
        names=[n.name for n in asset.names],
        description=asset.description or (asset_class.description if asset_class else None),
        seek_keys={sk.name: sk.value for sk in asset.seek_keys},
        parents=[AssetRow.of(p) for p in asset.parents],
        children=[AssetRow.of(c) for c in asset.children],
        archive_url=archive_url,
        files=files,
        pull_command=None if local else f"refgenie pull {target}",
        seek_command=f"refgenie seek {target}" if local else None,
    )


# --- Database layer tools ---


@_tool
def list_genomes(
    has_assets: bool | None = None,
    taxon_id: int | None = None,
    page_size: int | None = None,
    token: str | None = None,
) -> ResultPage[GenomeRow]:
    """List genomes, most assets first, then by species.

    Each row is small: digest, aliases, species and common name, asset count.
    20 rows per page by default, at most 100. Next: list_assets(alias=...) for
    a genome's files, or search_genomes(query) to find one by name.

    Args:
        has_assets: True for only genomes with assets (most have none).
        taxon_id: Only genomes of this NCBI taxon, e.g. 9606 for human.
        page_size: Rows per page (default 20, at most 100).
        token: The pagination.next_page_token of the previous page.
    """
    size = paging.clamp(page_size)
    fh = paging.filters_hash("list_genomes", has_assets=has_assets, taxon_id=taxon_id)
    page = _get_refgenie().genome.page(
        has_assets=has_assets,
        taxon_id=taxon_id,
        offset=paging.decode(token, fh),
        limit=size,
    )
    hint = _TAXON_HINT.format(taxon_id=taxon_id) if taxon_id and page.total == 0 else None
    return paging.result_page(page, [GenomeRow.of(g) for g in page.items], fh, hint)


@_tool
def search_genomes(
    query: str,
    has_assets: bool | None = None,
    taxon_id: int | None = None,
    page_size: int | None = None,
    token: str | None = None,
) -> ResultPage[GenomeRow]:
    """Search genomes by any name of their species (scientific, common, or another
    common name such as brown rat), alias, other genome names (GRCh38), digest,
    description, assembly source or assembly accession (case-insensitive
    substring match).

    Rows and paging are as in list_genomes (20 per page, at most 100). Matching
    is on names, not meaning. Next: list_assets(alias=...) for a genome's files.

    Args:
        query: The text to look for.
        has_assets: True for only genomes with assets (most have none).
        taxon_id: Only genomes of this NCBI taxon, e.g. 4577 for Zea mays.
        page_size: Rows per page (default 20, at most 100).
        token: The pagination.next_page_token of the previous page.
    """
    size = paging.clamp(page_size)
    fh = paging.filters_hash(
        "search_genomes", query=query, has_assets=has_assets, taxon_id=taxon_id
    )
    page = _get_refgenie().genome.page(
        search=Search(term=query),
        has_assets=has_assets,
        taxon_id=taxon_id,
        offset=paging.decode(token, fh),
        limit=size,
    )
    hint = None
    if page.total == 0:
        hint = (
            _TAXON_HINT.format(taxon_id=taxon_id)
            if taxon_id
            else _NO_MATCH_HINT.format(query=quote(query))
        )
    return paging.result_page(page, [GenomeRow.of(g) for g in page.items], fh, hint)


@_tool
def get_genome(alias: str | None = None, digest: str | None = None) -> GenomeDetail:
    """Get one genome's full record, including its FHR metadata and its species'
    names (scientific, common, other common names). Give exactly one of alias or
    digest. Next: list_assets with the same alias or digest.

    Args:
        alias: Genome alias, a name like 'hg38'.
        digest: Genome digest, 32 characters.
    """
    detail = _get_refgenie().genome.detail(_genome_digest(alias, digest))
    genome = detail.genome
    return GenomeDetail(
        **GenomeRow.of(detail).model_dump(),
        description=genome.description,
        synonyms=genome.synonyms,
        taxon_id=genome.taxon_id,
        assembly_source=genome.assembly_source,
        assembly_accession=genome.assembly_accession,
        assembly_level=genome.assembly_level,
        fhr=detail.fhr,
        species=SpeciesNames.of(detail.species),
    )


@_tool
def list_assets(
    alias: str | None = None,
    digest: str | None = None,
    asset_class: str | None = None,
    page_size: int | None = None,
    token: str | None = None,
) -> ResultPage[AssetRow]:
    """List the files you can get: assets, optionally of one genome and/or one asset class.

    Name the genome by alias or by digest, not both. Each row is small: digest,
    asset class, name, genome digest, size, serving modes. 20 rows per page by
    default, at most 100. Next: get_asset(digest) for files and the pull command.

    Args:
        alias: Optional genome alias (a name like 'hg38') to filter by.
        digest: Optional genome digest (32 characters) to filter by.
        asset_class: Optional asset group name to filter by (conventionally
            the asset class name, e.g. 'fasta' or 'bowtie2_index').
        page_size: Rows per page (default 20, at most 100).
        token: The pagination.next_page_token of the previous page.
    """
    size = paging.clamp(page_size)
    genome_digest = None if alias is None and digest is None else _genome_digest(alias, digest)
    fh = paging.filters_hash("list_assets", genome=genome_digest, asset_class=asset_class)
    page = _get_refgenie().asset.page(
        genome_digest=genome_digest,
        asset_group_name=asset_class,
        offset=paging.decode(token, fh),
        limit=size,
    )
    hint = None
    if genome_digest is not None and page.total == 0:
        hint = _NO_ASSETS_HINT.format(
            label=alias or digest, alias=quote(alias or ""), digest=genome_digest
        )
    return paging.result_page(page, [AssetRow.of(a) for a in page.items], fh, hint)


@_tool
def get_asset(digest: str) -> AssetDetail:
    """Get an asset's files and how to pull it, by its digest: file paths or URLs,
    the archive URL, the pull command, seek keys (relative to the asset folder),
    parents and children.

    Args:
        digest: The asset digest.
    """
    return _asset_detail(_get_refgenie().asset.get_by_digest(digest, with_links=True))


@_tool
def list_recipes(
    output_asset_class: str | None = None,
    page_size: int | None = None,
    token: str | None = None,
) -> ResultPage[RecipeRow]:
    """List recipes: how each asset class is built.

    Each row is small: name, version, output asset class, a one-line summary and
    what it needs. 20 rows per page by default, at most 100. Next:
    get_recipe(name) for the inputs and commands.

    Args:
        output_asset_class: Only recipes that build this asset class.
        page_size: Rows per page (default 20, at most 100).
        token: The pagination.next_page_token of the previous page.
    """
    size = paging.clamp(page_size)
    fh = paging.filters_hash("list_recipes", output_asset_class=output_asset_class)
    page = _get_refgenie().recipe.page(
        output_asset_class=output_asset_class, offset=paging.decode(token, fh), limit=size
    )
    return paging.result_page(page, [RecipeRow.of(r) for r in page.items], fh)


@_tool
def get_recipe(name: str, version: str | None = None) -> RecipeDetail:
    """Get one recipe in full: inputs, command templates, docker image, custom
    seek keys, output asset class, and the build command to run next.

    Args:
        name: The recipe name, from list_recipes.
        version: The recipe version; the latest if omitted.
    """
    return RecipeDetail.of(_get_refgenie().recipe.get(name, version))


@_tool
def list_asset_classes(
    page_size: int | None = None, token: str | None = None
) -> ResultPage[AssetClassRow]:
    """List asset classes: the kinds of asset (fasta, bowtie2_index, ...).

    Each row is small: name, version, a one-line summary and its seek key names.
    20 rows per page by default, at most 100. Next: list_recipes(output_asset_class=...)
    to see how one is built.

    Args:
        page_size: Rows per page (default 20, at most 100).
        token: The pagination.next_page_token of the previous page.
    """
    size = paging.clamp(page_size)
    fh = paging.filters_hash("list_asset_classes")
    page = _get_refgenie().asset_class.page(offset=paging.decode(token, fh), limit=size)
    return paging.result_page(page, [AssetClassRow.of(c) for c in page.items], fh)


@_tool
def lookup_digest(digest: str) -> DigestLookup:
    """Say what a digest is: a genome, an asset, or nothing here.

    Args:
        digest: A genome digest (32 characters) or an asset digest.
    """
    rgc = _get_refgenie()
    try:
        return DigestLookup(
            type="genome", digest=digest, genome=GenomeRow.of(rgc.genome.detail(digest))
        )
    except MissingGenomeError:
        pass
    try:
        return DigestLookup(
            type="asset", digest=digest, asset=AssetRow.of(rgc.asset.get_by_digest(digest))
        )
    except MissingAssetError:
        return DigestLookup(type="not_found", digest=digest)


@_tool
def catalog_summary() -> CatalogSummary:
    """See what this server holds: which server it is, counts of genomes, asset
    groups and assets, and the species that have assets. A good first call.
    Next: list_genomes(has_assets=true) or search_genomes(query)."""
    database = _get_refgenie().database
    species = [
        SpeciesCount(species=name, **counts) for name, counts in database.species_summary().items()
    ]
    species.sort(key=lambda s: s.assets, reverse=True)
    return CatalogSummary(server=_server_line(), **database.summary(), species_with_assets=species)


# --- RefgetStore layer tools ---


@_tool
def get_genome_sequences(
    alias: str | None = None,
    digest: str | None = None,
    name_contains: str | None = None,
    page_size: int | None = None,
    token: str | None = None,
) -> SequencePage:
    """List a genome's sequences, its chromosome sizes: name, length and refget
    digest, a page at a time. For a .chrom.sizes file, get_asset on the genome's
    fasta_index asset gives a download URL.

    Args:
        alias: Genome alias, a name like 'hg38'.
        digest: Genome digest, 32 characters.
        name_contains: Only sequences whose name contains this text (any case).
        page_size: Rows per page (default 100, at most 100).
        token: The pagination.next_page_token of the previous page.
    """
    genome = _genome_digest(alias, digest)
    level2 = _get_refgenie().genome.level2(genome)
    needle = (name_contains or "").lower()
    rows = [
        SequenceRow(name=n, length=length, digest=d)
        for n, length, d in zip(level2["names"], level2["lengths"], level2["sequences"])
        if needle in n.lower()
    ]
    size = paging.clamp(page_size or 100)
    fh = paging.filters_hash("get_genome_sequences", genome=genome, name_contains=name_contains)
    offset = paging.decode(token, fh)
    page = Page(items=rows[offset : offset + size], offset=offset, limit=size, total=len(rows))
    result = paging.result_page(page, page.items, fh)
    return SequencePage(
        **result.model_dump(),
        n_sequences=len(level2["names"]),
        total_length=sum(level2["lengths"]),
    )


@_tool
def get_sequence(
    name: str, start: int, end: int, alias: str | None = None, digest: str | None = None
) -> str:
    """Get the DNA bases of one region of a genome, without downloading the genome.

    Coordinates are 0-based and half-open, like BED: start=0, end=10 is the first
    10 bases. For 1-based coordinates from a genome browser (chr17:43,044,295-43,125,483),
    subtract 1 from the start only: start=43044294, end=43125483.
    At most 100,000 bases per call.
    name is the sequence name as this genome spells it (chr17 or 17);
    get_genome_sequences lists them.

    Args:
        name: Sequence name, as the genome spells it.
        start: First base, 0-based.
        end: End of the region, exclusive.
        alias: Genome alias, a name like 'hg38'.
        digest: Genome digest, 32 characters.
    """
    if end - start > MAX_SEQUENCE_REGION:
        raise ValueError(
            f"At most {MAX_SEQUENCE_REGION:,} bases per call; for longer regions run "
            "refgenie getseq or pull the genome's fasta asset."
        )
    return _get_refgenie().sequence.get(_genome_digest(alias, digest), name, start, end)


@_tool
def compare_genomes(
    alias_a: str | None = None,
    digest_a: str | None = None,
    alias_b: str | None = None,
    digest_b: str | None = None,
) -> GenomeComparison:
    """Say whether two genomes are the same and how they differ, in one sentence,
    with a map of renamed sequences (chrI->I).

    comparison holds the GA4GH seqcol comparison object: digests, attributes
    (a_only, b_only, a_and_b) and array_elements (a_count, b_count,
    a_and_b_count, a_and_b_same_order, each keyed by attribute).

    Name genome A with exactly one of alias_a or digest_a, and genome B with
    exactly one of alias_b or digest_b.

    Args:
        alias_a: First genome's alias, a name like 'hg38'.
        digest_a: First genome's digest, 32 characters.
        alias_b: Second genome's alias.
        digest_b: Second genome's digest.
    """
    a = _genome_digest(alias_a, digest_a, "_a")
    b = _genome_digest(alias_b, digest_b, "_b")
    side_a, side_b = _get_refgenie().genome.compare_sides(a, b)
    comparison = compare_level2(side_a, side_b, a, b)
    summary, name_map = summarize_comparison(side_a, side_b, comparison, alias_a or a, alias_b or b)
    return GenomeComparison(summary=summary, name_map=name_map, comparison=comparison)
