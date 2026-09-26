"""
MCP tool definitions for Refgenie.

All tools are read-only queries against the Refgenie object.
Shared between stdio (local) and Streamable HTTP (server) transports.
"""

import json

try:
    # mcp 2.x moved the high-level server to mcp.server.mcpserver.MCPServer
    # (was mcp.server.fastmcp.FastMCP).
    from mcp.server.mcpserver import MCPServer
except ImportError as e:
    raise ImportError(
        "The 'mcp' extras are not installed. Please install refgenie with the 'mcp' extras to use the MCP server."
    ) from e

from sqlmodel import Session, or_, select

from refgenie import __version__ as _refgenie_version
from refgenie.db.tables import Genome
from refgenie.exceptions import MissingAssetError
from refgenie.models import GenomeAlias, GenomeDigest
from refgenie.utils.search import SearchOperator, alias_digests_matching

mcp = MCPServer("refgenie", version=_refgenie_version)

_rgc = None


def _get_refgenie():
    """Lazy singleton. Overridden by set_refgenie() for server mode."""
    global _rgc
    if _rgc is None:
        from refgenie import Refgenie

        _rgc = Refgenie()
    return _rgc


def set_refgenie(rgc):
    """Inject an existing Refgenie instance (used by server/main.py)."""
    global _rgc
    _rgc = rgc


def _json_response(data) -> str:
    """Serialize response data to JSON string."""
    return json.dumps(data, default=str)


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
        return GenomeDigest(digest)
    return _get_refgenie().alias.resolve(GenomeAlias(alias))


def _genome_aliases(genome_digest: GenomeDigest) -> list[str]:
    """Alias names for a genome, from the mode-selected alias manager.

    Reads from rgc.alias (store-backed in local mode, SQL-backed in server
    mode) rather than the SQL Genome.aliases relationship, which is empty in
    local mode. Mirrors the server router (server/routers/catalog.py).
    """
    return _get_refgenie().alias.get_for_genome(genome_digest)


def _genome_summary(genome: Genome) -> dict:
    """The genome payload shared by list_genomes and search_genomes."""
    return {
        "digest": genome.digest,
        "aliases": _genome_aliases(genome.digest),
        "species_name": genome.species_name,
        "common_name": genome.common_name,
        "taxon_id": genome.taxon_id,
        "assembly_source": genome.assembly_source,
        "assembly_accession": genome.assembly_accession,
        "description": genome.description,
    }


# --- Database layer tools ---


@mcp.tool()
def list_genomes() -> str:
    """List all genomes with their aliases, species, and description."""
    rgc = _get_refgenie()
    return _json_response([_genome_summary(g) for g in rgc.genome.list_all()])


@mcp.tool()
def list_asset_classes() -> str:
    """List all registered asset classes with their serving modes."""
    rgc = _get_refgenie()
    asset_classes = rgc.asset_class.list_all()
    result = []
    for ac in asset_classes:
        result.append(
            {
                "name": ac.name,
                "version": ac.version,
                "description": ac.description,
                "serving_modes": ac.serving_modes,
            }
        )
    return _json_response(result)


@mcp.tool()
def list_recipes() -> str:
    """List all registered recipes with their inputs and command templates."""
    rgc = _get_refgenie()
    recipes = rgc.recipe.list_all()
    result = []
    for r in recipes:
        result.append(
            {
                "name": r.name,
                "version": r.version,
                "description": r.description,
                "input_params": r.input_params,
                "input_files": r.input_files,
                "input_assets": r.input_assets,
                "command_templates": r.command_templates,
                "docker_image": r.docker_image,
            }
        )
    return _json_response(result)


@mcp.tool()
def list_assets(
    alias: str | None = None, digest: str | None = None, asset_class: str | None = None
) -> str:
    """List assets, optionally filtered by one genome and/or asset class.

    Name the genome by alias or by digest, not both.

    Args:
        alias: Optional genome alias (a name like 'hg38') to filter by.
        digest: Optional genome digest (32 characters) to filter by.
        asset_class: Optional asset group name to filter by (conventionally
            the asset class name).
    """
    rgc = _get_refgenie()
    genome_digests = None if alias is None and digest is None else [_genome_digest(alias, digest)]
    assets = rgc.asset.list_assets(
        genome_digests=genome_digests,
        asset_group_name=asset_class,
    )
    result = []
    for a in assets:
        result.append(
            {
                "digest": a.digest,
                "name": a.name,
                "asset_group": a.asset_group.name if a.asset_group else None,
                "genome": a.asset_group.genome_digest if a.asset_group else None,
                "size": a.size,
                "seek_keys": {sk.name: sk.value for sk in a.seek_keys} if a.seek_keys else {},
            }
        )
    return _json_response(result)


@mcp.tool()
def search_genomes(query: str) -> str:
    """Search genomes by species/common name, alias, description, assembly source,
    or assembly accession using substring matching.

    Args:
        query: Search string matched against genome species/common names, aliases,
            descriptions, assembly source, and assembly accession.
    """
    rgc = _get_refgenie()
    # Matched in SQL, not in Python: this tool runs against the full hosted
    # catalog, which must never be pulled into memory to answer a search.
    # Alias matching goes through the mode-selected alias manager instead of a
    # SQL join, since aliases may live in the store rather than the SQL table.
    pattern = f"%{query}%"
    alias_matches = alias_digests_matching(rgc.alias, query, SearchOperator.CONTAINS)
    statement = select(Genome).where(
        or_(
            Genome.species_name.ilike(pattern),
            Genome.common_name.ilike(pattern),
            Genome.description.ilike(pattern),
            Genome.assembly_source.ilike(pattern),
            Genome.assembly_accession.ilike(pattern),
            Genome.digest.in_(alias_matches),
        )
    )
    with Session(rgc.database_engine, expire_on_commit=False) as session:
        genomes = session.exec(statement).unique().all()
    return _json_response([_genome_summary(g) for g in genomes])


@mcp.tool()
def get_genome(alias: str | None = None, digest: str | None = None) -> str:
    """Get detailed info for a single genome. Give exactly one of alias or digest.

    Args:
        alias: Genome alias, a name like 'hg38'.
        digest: Genome digest, 32 characters.
    """
    rgc = _get_refgenie()
    digest = _genome_digest(alias, digest)
    genome = rgc.genome.get(digest)
    asset_groups = []
    for ag in genome.asset_groups:
        assets = []
        for a in ag.assets:
            assets.append(
                {
                    "digest": a.digest,
                    "name": a.name,
                    "size": a.size,
                }
            )
        asset_groups.append(
            {
                "name": ag.name,
                "description": ag.description,
                "assets": assets,
            }
        )

    result = {
        "digest": genome.digest,
        "aliases": _genome_aliases(digest),
        "species_name": genome.species_name,
        "common_name": genome.common_name,
        "taxon_id": genome.taxon_id,
        "assembly_source": genome.assembly_source,
        "assembly_accession": genome.assembly_accession,
        "assembly_level": genome.assembly_level,
        "description": genome.description,
        "asset_groups": asset_groups,
    }
    return _json_response(result)


@mcp.tool()
def get_asset(digest: str) -> str:
    """Get detailed info for a single asset by its digest.

    Args:
        digest: The asset digest.
    """
    rgc = _get_refgenie()
    asset = rgc.asset.get_by_digest(digest)
    if asset is None:
        raise MissingAssetError(digest=digest)
    result = {
        "digest": asset.digest,
        "name": asset.name,
        "description": asset.description,
        "size": asset.size,
        "path": asset.path,
        "asset_group": asset.asset_group.name if asset.asset_group else None,
        "genome": asset.asset_group.genome_digest if asset.asset_group else None,
        "seek_keys": {sk.name: sk.value for sk in asset.seek_keys} if asset.seek_keys else {},
        "parents": [p.digest for p in asset.parents] if asset.parents else [],
        "children": [c.digest for c in asset.children] if asset.children else [],
    }
    return _json_response(result)


@mcp.tool()
def lookup_digest(digest: str) -> str:
    """Universal digest lookup - tries genome digest, then asset digest.

    Args:
        digest: A digest to look up (could be genome or asset).
    """
    rgc = _get_refgenie()

    if rgc.genome.exists(digest):
        genome = rgc.genome.get(digest)
        return _json_response(
            {
                "type": "genome",
                "digest": genome.digest,
                "aliases": _genome_aliases(genome.digest),
                "species_name": genome.species_name,
                "description": genome.description,
            }
        )

    asset = rgc.asset.get_by_digest(digest)
    if asset is None:
        return _json_response({"type": "not_found", "digest": digest})
    return _json_response(
        {
            "type": "asset",
            "digest": asset.digest,
            "name": asset.name,
            "asset_group": asset.asset_group.name if asset.asset_group else None,
            "genome": asset.asset_group.genome_digest if asset.asset_group else None,
        }
    )


# --- RefgetStore layer tools ---


@mcp.tool()
def get_genome_metadata(alias: str | None = None, digest: str | None = None) -> str:
    """Get sequence collection metadata: number of sequences, total length, source.

    Give exactly one of alias or digest.

    Args:
        alias: Genome alias, a name like 'hg38'.
        digest: Genome digest, 32 characters.
    """
    rgc = _get_refgenie()
    result = rgc.genome.get_metadata(_genome_digest(alias, digest))
    return _json_response(result)


@mcp.tool()
def get_genome_sequences(alias: str | None = None, digest: str | None = None) -> str:
    """Get sequence collection level2 data: names, lengths, and sequence digests.

    Give exactly one of alias or digest.

    Args:
        alias: Genome alias, a name like 'hg38'.
        digest: Genome digest, 32 characters.
    """
    rgc = _get_refgenie()
    level2 = rgc.store_router.get_collection_level2(_genome_digest(alias, digest))
    return _json_response(level2)


@mcp.tool()
def compare_genomes(
    alias_a: str | None = None,
    digest_a: str | None = None,
    alias_b: str | None = None,
    digest_b: str | None = None,
) -> str:
    """Compare two genomes using seqcol comparison.

    Name genome A with exactly one of alias_a or digest_a, and genome B with
    exactly one of alias_b or digest_b.

    Args:
        alias_a: First genome's alias, a name like 'hg38'.
        digest_a: First genome's digest, 32 characters.
        alias_b: Second genome's alias.
        digest_b: Second genome's digest.
    """
    rgc = _get_refgenie()
    result = rgc.genome.compare(
        _genome_digest(alias_a, digest_a, "_a"), _genome_digest(alias_b, digest_b, "_b")
    )
    return _json_response(result)
