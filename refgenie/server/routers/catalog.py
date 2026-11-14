from dataclasses import asdict
from typing import Any

from fastapi import APIRouter, Depends, HTTPException, Query

from refgenie.core import Refgenie
from refgenie.db.tables import (
    AliasPublic,
    AssetClassPublic,
    AssetGroupPublic,
    ConfigurationPublic,
    StagedAssetPublic,
)
from refgenie.exceptions import (
    MissingAliasError,
    MissingAssetError,
    MissingConfigDataError,
    MissingGenomeError,
    UnservedStoreError,
)
from refgenie.models import GenomeAlias, display_alias
from refgenie.server.dependencies import DbId, get_refgenie
from refgenie.server.errors import ErrorCode
from refgenie.server.pagination import PageParams, PageResponse, get_page_params, respond
from refgenie.server.routers import helpers
from refgenie.server.schemas import (
    AliasResponse,
    AssetClassResponse,
    AssetResponse,
    GenomeDetailResponse,
    GenomeResponse,
    GenomeSequencesResponse,
    RecipeResponse,
    SpeciesResponse,
)
from refgenie.server.search import SearchParams, get_search_params, invalid_search_is_422


def _taxon_uri(taxon_id: int | None) -> str | None:
    """identifiers.org taxonomy URI derived from an NCBI taxon id, or None."""
    if taxon_id is None:
        return None
    return f"https://identifiers.org/taxonomy:{taxon_id}"


def _asset_response(asset) -> AssetResponse:
    return AssetResponse.model_validate(asset, from_attributes=True)


router = APIRouter()


@router.get("/genomes", response_model=PageResponse[GenomeResponse])
def list_genomes(
    digest: str | None = Query(None, description="The genome digest"),
    alias: str | None = Query(None, description="The genome alias"),
    has_assets: bool | None = Query(
        None, description="Only genomes with (true) or without (false) assets"
    ),
    taxon_id: int | None = Query(None, description="Only genomes of this NCBI taxon id"),
    store: str | None = Query(
        None, description="Only genomes present in these stores (comma-separated store names)"
    ),
    search: SearchParams = Depends(get_search_params),
    paging: PageParams = Depends(get_page_params),
    rgc: Refgenie = Depends(get_refgenie),
):
    """
    List genomes with their aliases and asset counts, most asset groups first,
    then by species name and digest. With ``q``, best matches come first:
    an exact name, then a whole word of a name, then a word start, then the rest.

    Exact-match filters:
    - digest: Filter by exact genome digest
    - alias: Filter by exact genome alias
    - has_assets: Only genomes with (true) or without (false) assets
    - taxon_id: Filter by NCBI taxon id
    - store: Only genomes present in any of these stores (comma-separated).
      This is membership, so a genome held by several stores matches each.
      An unknown store name is a 422.

    Searchable fields: see ``SEARCH_FIELDS`` in ``refgenie.db.queries.genome``.
    An unknown alias gives an empty page, not a 404.
    """
    stores = [s.strip() for s in (store or "").split(",") if s.strip()] or None
    # UnknownStoreError is a ValueError, so this also turns it into a 422.
    with invalid_search_is_422():
        page = rgc.genome.page(
            digest=digest or None,
            alias=GenomeAlias(alias) if alias else None,
            search=search.to_search(),
            has_assets=has_assets,
            taxon_id=taxon_id,
            stores=stores,
            offset=paging.offset,
            limit=paging.page_size,
        )
    return respond(page, paging, GenomeResponse.of)


def _genome_not_found(digest: str, rgc: Refgenie) -> HTTPException:
    """The 404 for a digest with no genome row.

    In local mode, aliases live in the refget store and genomes in SQL, so an
    interrupted ``genome remove`` can leave aliases pointing at a missing
    genome; that case gets a ``stale_alias`` error naming the fix.
    """
    try:
        stale = rgc.alias.get_for_genome(digest)
    except MissingGenomeError:
        stale = []
    if not stale:
        return HTTPException(status_code=404, detail=f"Genome {digest} not found")
    names = ", ".join(sorted(stale))
    first = sorted(stale)[0]
    return HTTPException(
        status_code=404,
        detail={
            "code": ErrorCode.STALE_ALIAS,
            "message": (
                f"No genome record for {digest}, but the alias(es) {names} still point at it. "
                "The refget store and the catalog disagree, which a genome removal interrupted "
                f"partway through will do. Remove the stale alias with `refgenie alias remove {first}`."
            ),
        },
    )


def _unserved(error: UnservedStoreError) -> HTTPException:
    """The 404 for a genome whose store this server does not serve."""
    return HTTPException(
        status_code=404, detail={"code": ErrorCode.GENOME_NOT_FOUND, "message": str(error)}
    )


def _no_collection(error: Exception, digests: str) -> HTTPException:
    """The 404 for a genome with no sequence collection in its store."""
    if isinstance(error, UnservedStoreError):
        return _unserved(error)
    return HTTPException(
        status_code=404,
        detail={
            "code": ErrorCode.GENOME_NOT_FOUND,
            "message": f"No sequence collection for genome {digests} in this refgenie's stores",
        },
    )


@router.get("/genomes/{digest}", response_model=GenomeDetailResponse)
def get_genome(digest: str, rgc: Refgenie = Depends(get_refgenie)):
    """Return a genome's core columns plus its FHR metadata, taxon URI and species."""
    try:
        detail = rgc.genome.detail(digest)
    except UnservedStoreError as e:
        raise _unserved(e) from None
    except MissingGenomeError:
        raise _genome_not_found(digest, rgc) from None
    response = GenomeDetailResponse.model_validate(detail.genome, from_attributes=True)
    response.common_name = detail.common_name
    response.taxon_uri = _taxon_uri(detail.genome.taxon_id)
    response.fhr = detail.fhr
    response.species = SpeciesResponse.of(detail.species)
    response.refget_store_url = detail.store_url
    return response


@router.get("/genomes/{digest}/sequences", response_model=GenomeSequencesResponse)
def get_genome_sequences(digest: str, rgc: Refgenie = Depends(get_refgenie)):
    """Return a genome's sequences (names, lengths, digests, aliases) from the refget store that holds it."""
    try:
        collection = rgc.genome.sequences(digest)
    except MissingGenomeError as e:
        raise _no_collection(e, digest) from None
    return GenomeSequencesResponse.model_validate(
        {
            "digest": collection.digest,
            "n_sequences": len(collection.sequences),
            "total_length": sum(row.length for row in collection.sequences),
            "level1": collection.level1,
            "collection_aliases": [asdict(a) for a in collection.collection_aliases],
            "sequences": [asdict(row) for row in collection.sequences],
        }
    )


@router.get("/genomes/{digest}/compare/{other}", response_model=dict[str, Any])
def compare_genomes(digest: str, other: str, rgc: Refgenie = Depends(get_refgenie)):
    """Compare two genomes' sequence collections (GA4GH seqcol comparison).

    Adds ``summary``, one sentence on how they differ, and ``name_map``, A's
    name to B's for each sequence they share under another name, so a client
    can word the result without fetching both collections.
    """
    from refgenie.utils.seqcol import compare_level2, seqcol_level2, summarize_comparison

    try:
        a = rgc.genome.get_served(digest)
        b = rgc.genome.get_served(other)
        # The loaded stores, never a genome's remote_url; the same comparison
        # `refgenie compare` and MCP return.
        side_a = seqcol_level2(rgc.store_router.store_for_genome(digest, a.store_name), digest)
        side_b = seqcol_level2(rgc.store_router.store_for_genome(other, b.store_name), other)
    except (MissingGenomeError, ValueError) as e:  # ValueError: no level 2 in the store
        raise _no_collection(e, f"{digest} or {other}") from None
    result = compare_level2(side_a, side_b, digest, other)
    aliases = rgc.alias.get_for_genomes([digest, other])
    label_a, label_b = (display_alias(aliases.get(d, [])) or d for d in (digest, other))
    summary, name_map = summarize_comparison(side_a, side_b, result, label_a, label_b)
    return {**result, "summary": summary, "name_map": name_map}


@router.get("/asset_groups", response_model=PageResponse[AssetGroupPublic])
def list_asset_groups(
    genome_digest: str | None = Query(None, description="The genome digest"),
    asset_class: str | None = Query(None, description="The asset class"),
    asset_group_name: str | None = Query(None, description="The asset group name"),
    asset_group_id: DbId | None = Query(None, description="The asset group ID"),
    search: SearchParams = Depends(get_search_params),
    paging: PageParams = Depends(get_page_params),
    rgc: Refgenie = Depends(get_refgenie),
):
    """
    List asset groups, by genome digest, name and id.

    Exact-match filters:
    - genome_digest: Filter by exact genome digest
    - asset_class: Filter by exact asset class name
    - asset_group_name: Filter by exact asset group name
    - asset_group_id: Filter by exact asset group ID

    Searchable columns: name
    """
    with invalid_search_is_422():
        page = rgc.asset.group.page(
            genome_digest=genome_digest,
            asset_class=asset_class,
            name=asset_group_name,
            id=asset_group_id,
            search=search.to_search(),
            offset=paging.offset,
            limit=paging.page_size,
        )
    return respond(page, paging, lambda g: AssetGroupPublic.model_validate(g, from_attributes=True))


@router.get("/asset_groups/{id}", response_model=AssetGroupPublic)
def get_asset_group(id: DbId, rgc: Refgenie = Depends(get_refgenie)):
    asset_group = rgc.asset.group.get_by_id(id)
    if not asset_group:
        raise HTTPException(status_code=404, detail=f"Asset group {id} not found")
    return asset_group


@router.get("/assets/{digest}", response_model=AssetResponse)
def get_asset(digest: str, rgc: Refgenie = Depends(get_refgenie)):
    try:
        return _asset_response(rgc.asset.get_by_digest(digest))
    except MissingAssetError:
        raise HTTPException(status_code=404, detail=f"Asset {digest} not found") from None


@router.get("/assets", response_model=PageResponse[AssetResponse])
def list_assets(
    name: str | None = Query(None, description="The tag name"),
    asset_group_name: str | None = Query(None, description="The asset group name"),
    genome_digest: str | None = Query(None, description="The genome id"),
    recipe_name: str | None = Query(
        None, description="The recipe name that was used to create the asset"
    ),
    asset_group_id: DbId | None = Query(None, description="The asset group ID to filter assets by"),
    search: SearchParams = Depends(get_search_params),
    paging: PageParams = Depends(get_page_params),
    rgc: Refgenie = Depends(get_refgenie),
):
    """
    List assets, by genome digest, asset group name and digest.

    Exact-match filters:
    - name: Filter by exact asset name (any of the asset's names)
    - asset_group_name: Filter by exact asset group name
    - genome_digest: Filter by exact genome digest
    - recipe_name: Filter by exact recipe name
    - asset_group_id: Filter by exact asset group ID

    Searchable columns: name, digest, path
    """
    with invalid_search_is_422():
        page = rgc.asset.page(
            genome_digest=genome_digest,
            asset_group_name=asset_group_name,
            asset_group_id=asset_group_id,
            name=name,
            recipe_name=recipe_name,
            search=search.to_search(),
            offset=paging.offset,
            limit=paging.page_size,
        )
    return respond(page, paging, _asset_response)


@router.get("/assets/{asset_digest}/files")
def list_asset_files(asset_digest: str, rgc: Refgenie = Depends(get_refgenie)):
    """
    List files available for individual download from an asset.

    Returns ``{"asset_digest": ..., "files": [...]}`` from the asset's
    file-mode staged record.
    """
    _, sa = helpers.require_published(rgc, asset_digest, "file")

    return {"asset_digest": asset_digest, "files": sa.directory_contents}


@router.get("/asset_classes", response_model=PageResponse[AssetClassPublic])
def list_asset_classes(
    name: str | None = Query(None, description="The asset class name"),
    version: str | None = Query(None, description="The asset class version"),
    search: SearchParams = Depends(get_search_params),
    paging: PageParams = Depends(get_page_params),
    rgc: Refgenie = Depends(get_refgenie),
):
    """
    List asset classes, by name and version.

    Exact-match filters:
    - name: Filter by exact asset class name
    - version: Filter by exact asset class version

    Searchable columns: name, version, description
    """
    with invalid_search_is_422():
        page = rgc.asset_class.page(
            name=name,
            version=version,
            search=search.to_search(),
            offset=paging.offset,
            limit=paging.page_size,
        )
    return respond(page, paging, lambda c: AssetClassPublic.model_validate(c, from_attributes=True))


@router.get("/asset_classes/{id}", response_model=AssetClassResponse)
def get_asset_class(id: DbId, rgc: Refgenie = Depends(get_refgenie)):
    """One asset class, with its seek keys."""
    asset_class = rgc.asset_class.get_by_id(id)
    if not asset_class:
        raise HTTPException(status_code=404, detail=f"Asset class {id} not found")
    return AssetClassResponse.model_validate(asset_class, from_attributes=True)


@router.get("/recipes", response_model=PageResponse[RecipeResponse])
def list_recipes(
    output_asset_class: str | None = Query(None, description="The asset class"),
    name: str | None = Query(None, description="The recipe name"),
    version: str | None = Query(None, description="The recipe version"),
    search: SearchParams = Depends(get_search_params),
    paging: PageParams = Depends(get_page_params),
    rgc: Refgenie = Depends(get_refgenie),
):
    """
    List recipes, by name and version.

    Exact-match filters:
    - output_asset_class: Filter by exact output asset class name
    - name: Filter by exact recipe name
    - version: Filter by exact recipe version

    Searchable columns: name, version, description
    """
    with invalid_search_is_422():
        page = rgc.recipe.page(
            output_asset_class=output_asset_class,
            name=name,
            version=version,
            search=search.to_search(),
            offset=paging.offset,
            limit=paging.page_size,
        )
    return respond(page, paging, RecipeResponse.of)


@router.get("/recipes/{id}", response_model=RecipeResponse)
def get_recipe(id: DbId, rgc: Refgenie = Depends(get_refgenie)):
    recipe = rgc.recipe.get_by_id(id)
    if not recipe:
        raise HTTPException(status_code=404, detail=f"Recipe {id} not found")
    return RecipeResponse.of(recipe)


@router.get("/configurations", response_model=PageResponse[ConfigurationPublic])
def list_configurations(
    paging: PageParams = Depends(get_page_params),
    rgc: Refgenie = Depends(get_refgenie),
):
    page = rgc.configuration.page(offset=paging.offset, limit=paging.page_size)
    return respond(
        page, paging, lambda c: ConfigurationPublic.model_validate(c, from_attributes=True)
    )


@router.get("/configurations/{id}", response_model=ConfigurationPublic)
def get_configuration(id: DbId, rgc: Refgenie = Depends(get_refgenie)):
    try:
        return rgc.configuration.get(id)
    except MissingConfigDataError:
        raise HTTPException(status_code=404, detail=f"Configuration {id} not found") from None


@router.get("/staged_assets", response_model=PageResponse[StagedAssetPublic])
def list_staged_assets(
    asset_digest: str | None = Query(None, description="The asset digest"),
    mode: str | None = Query(None, description="The staging mode (file or archive)"),
    search: SearchParams = Depends(get_search_params),
    paging: PageParams = Depends(get_page_params),
    rgc: Refgenie = Depends(get_refgenie),
):
    """
    List staged assets, by asset digest and mode.

    Searchable columns: asset_digest, mode
    """
    with invalid_search_is_422():
        page = rgc.stage.page(
            asset_digest=asset_digest,
            mode=mode,
            search=search.to_search(),
            offset=paging.offset,
            limit=paging.page_size,
        )
    return respond(
        page, paging, lambda s: StagedAssetPublic.model_validate(s, from_attributes=True)
    )


@router.get("/staged_assets/{id}", response_model=StagedAssetPublic)
def get_staged_asset(id: DbId, rgc: Refgenie = Depends(get_refgenie)):
    staged = rgc.stage.get(id)
    if not staged:
        raise HTTPException(status_code=404, detail=f"Staged asset {id} not found")
    return staged


@router.get("/relationships/{asset_digest}")
def get_asset_relationships(
    asset_digest: str,
    expand: bool = Query(
        False,
        description="If true, return the full asset objects instead of bare digests",
    ),
    rgc: Refgenie = Depends(get_refgenie),
):
    """
    Get parent and child relationships for a specific asset.

    Returns a dictionary with 'parents' and 'children' keys containing
    either asset digests (expand=false) or full asset objects (expand=true).
    """
    try:
        asset = rgc.asset.get_by_digest(asset_digest, with_links=True)
    except MissingAssetError:
        raise HTTPException(status_code=404, detail=f"Asset {asset_digest} not found") from None

    if expand:
        parents = [_asset_response(parent) for parent in asset.parents]
        children = [_asset_response(child) for child in asset.children]
    else:
        parents = [parent.digest for parent in asset.parents if parent.digest]
        children = [child.digest for child in asset.children if child.digest]

    return {"asset_digest": asset_digest, "parents": parents, "children": children}


@router.get("/aliases", response_model=PageResponse[AliasPublic])
def list_aliases(
    name: str | None = Query(None, description="The alias name"),
    genome_digest: str | None = Query(None, description="The genome digest"),
    search: SearchParams = Depends(get_search_params),
    paging: PageParams = Depends(get_page_params),
    rgc: Refgenie = Depends(get_refgenie),
):
    """
    List aliases, by name.

    Exact-match filters:
    - name: Filter by exact alias name
    - genome_digest: Filter by exact genome digest

    Searchable columns: name
    """
    with invalid_search_is_422():
        page = rgc.alias.page(
            name=name,
            genome_digest=genome_digest,
            search=search.to_search(),
            offset=paging.offset,
            limit=paging.page_size,
        )
    return respond(page, paging, lambda r: AliasPublic(name=r.name, genome_digest=r.genome_digest))


@router.get("/aliases/{name}", response_model=AliasResponse)
def get_alias(name: str, rgc: Refgenie = Depends(get_refgenie)):
    """Resolve a genome alias to its digest, sequence collection, and FHR metadata."""
    try:
        digest = rgc.alias.resolve(GenomeAlias(name))
    except MissingAliasError:
        raise HTTPException(status_code=404, detail=f"Alias '{name}' not found")

    try:
        genome = rgc.genome.get_served(digest)
        collection = rgc.store_router.get_collection_level2(digest, genome.store_name)
    except UnservedStoreError as e:
        raise _unserved(e) from None
    except MissingGenomeError:
        raise HTTPException(
            status_code=404, detail=f"Collection data not found for {digest}"
        ) from None

    _fhr = rgc.store_router.get_fhr_metadata(digest, genome.store_name)
    fhr = _fhr.to_dict() if _fhr is not None else None

    return AliasResponse(
        alias=name,
        digest=digest,
        source="server",
        collection=collection,
        fhr=fhr,
    )
