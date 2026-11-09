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
from pathlib import Path

try:
    from mcp.server.mcpserver import MCPServer
    from mcp.server.mcpserver.exceptions import ToolError
except ImportError as e:
    raise ImportError(
        "The 'mcp' extras are not installed. Please install refgenie with the 'mcp' extras to use the MCP server."
    ) from e

from refgenie import __version__ as _refgenie_version
from refgenie.db.tables import Asset
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
    GenomeDetail,
    GenomeRow,
    RecipeDetail,
    RecipeRow,
    SpeciesCount,
    SpeciesNames,
)
from refgenie.models import GenomeAlias, GenomeDigest
from refgenie.utils.search import Search

INSTRUCTIONS = """\
refgenie is a catalog of reference genomes and the assets built from them (FASTA files, aligner indexes, annotations).
A genome has a 32-character digest and any number of aliases (names like hg38). Tools take either alias= or digest=, never a guess between them.
Most genomes have no assets; pass has_assets=true to list only genomes with files.
The usual path: search_genomes or list_genomes, then list_assets(alias=...), then get_asset(digest) for files and the pull command.
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

_rgc = None
_serves_downloads = False

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
    "has_assets=true shows only genomes with files."
)


def _get_refgenie():
    """Lazy singleton. Overridden by set_refgenie() for server mode."""
    global _rgc
    if _rgc is None:
        from refgenie import Refgenie

        _rgc = Refgenie()
    return _rgc


def set_refgenie(rgc, *, serves_downloads: bool = False):
    """Inject an existing Refgenie instance (used by server/main.py).

    Args:
        rgc: The Refgenie instance the tools read.
        serves_downloads: Whether this process also serves the HTTP download
            routes, so ``get_asset`` gives download URLs instead of a local path.
    """
    global _rgc, _serves_downloads
    _rgc = rgc
    _serves_downloads = serves_downloads


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
    downloads = local_path = None
    if _serves_downloads:
        downloads = {}
        if "archive" in asset.serving_modes:
            downloads["archive"] = f"/v4/archives/{asset.digest}/download"
        if "file" in asset.serving_modes:
            downloads["file"] = f"/v4/assets/{asset.digest}/files"
    elif asset.path is not None:
        local_path = str(Path(rgc.genome_folder) / asset.path)
    return AssetDetail(
        **AssetRow.of(asset).model_dump(),
        asset_group=asset.asset_group_name,
        names=[n.name for n in asset.names],
        description=asset.description,
        seek_keys={sk.name: sk.value for sk in asset.seek_keys},
        parents=[AssetRow.of(p) for p in asset.parents],
        children=[AssetRow.of(c) for c in asset.children],
        downloads=downloads,
        local_path=local_path,
        pull_command=f"refgenie pull {aliases[0]}/{asset.asset_group_name}:{asset.name}"
        if aliases
        else None,
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
    return paging.result_page(page, [GenomeRow.of(g) for g in page.items], fh)


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
    hint = _NO_MATCH_HINT if page.total == 0 else None
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
    """List assets, optionally of one genome and/or one asset class.

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
    return paging.result_page(page, [AssetRow.of(a) for a in page.items], fh)


@_tool
def get_asset(digest: str) -> AssetDetail:
    """Get one asset's full record by its digest: names, seek keys, parents and
    children, and how to get it (download URLs on a server, the local path
    otherwise, and the refgenie pull command).

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
    """Count the catalog: genomes, asset groups and assets, and the species that
    have assets. A good first call to see what this server holds."""
    database = _get_refgenie().database
    species = [
        SpeciesCount(species=name, **counts) for name, counts in database.species_summary().items()
    ]
    species.sort(key=lambda s: s.assets, reverse=True)
    return CatalogSummary(**database.summary(), species_with_assets=species)


# --- RefgetStore layer tools ---


@_tool
def get_genome_metadata(alias: str | None = None, digest: str | None = None) -> dict:
    """Get sequence collection metadata: number of sequences, total length, source.

    Give exactly one of alias or digest.

    Args:
        alias: Genome alias, a name like 'hg38'.
        digest: Genome digest, 32 characters.
    """
    return _get_refgenie().genome.get_metadata(_genome_digest(alias, digest))


@_tool
def get_genome_sequences(alias: str | None = None, digest: str | None = None) -> dict:
    """Get sequence collection level2 data: names, lengths, and sequence digests.

    Give exactly one of alias or digest.

    Args:
        alias: Genome alias, a name like 'hg38'.
        digest: Genome digest, 32 characters.
    """
    return _get_refgenie().genome.level2(_genome_digest(alias, digest))


@_tool
def compare_genomes(
    alias_a: str | None = None,
    digest_a: str | None = None,
    alias_b: str | None = None,
    digest_b: str | None = None,
) -> dict:
    """Compare two genomes using seqcol comparison.

    Returns the GA4GH seqcol comparison object: digests, attributes (a_only,
    b_only, a_and_b) and array_elements (a_count, b_count, a_and_b_count,
    a_and_b_same_order, each keyed by attribute).

    Name genome A with exactly one of alias_a or digest_a, and genome B with
    exactly one of alias_b or digest_b.

    Args:
        alias_a: First genome's alias, a name like 'hg38'.
        digest_a: First genome's digest, 32 characters.
        alias_b: Second genome's alias.
        digest_b: Second genome's digest.
    """
    return _get_refgenie().genome.compare(
        _genome_digest(alias_a, digest_a, "_a"), _genome_digest(alias_b, digest_b, "_b")
    )
