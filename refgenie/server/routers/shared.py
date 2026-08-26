from fastapi import APIRouter, Depends, HTTPException, Query
from sqlalchemy import or_
from sqlalchemy.orm import selectinload
from sqlmodel import Session, select

from refgenie.server.dependencies import get_db_session, get_refgenie
from refgenie.server.routers import helpers
from refgenie.exceptions import MissingAliasError
from refgenie.core import Refgenie
from refgenie.db.tables import (
    StagedAsset,
    StagedAssetPublic,
    AliasPublic,
    Asset,
    AssetClass,
    AssetClassPublic,
    AssetGroup,
    AssetGroupPublic,
    AssetName,
    Configuration,
    ConfigurationPublic,
    Genome,
    Recipe,
    RecipePublic,
)
from refgenie.utils.pagination import (
    PaginatedResponse,
    PaginationMeta,
    PaginationParams,
    get_pagination_params,
    paginate_list,
    paginate_query,
)
from refgenie.utils.search import (
    SearchParams,
    _build_search_condition,
    _validate_search_fields,
    alias_digests_matching,
    apply_search_to_query,
    name_matches_operator,
    get_search_params,
)

from refgenie.server.schemas import (
    AliasResponse,
    AssetResponse,
    GenomeDetailResponse,
    GenomeResponse,
)


def _taxon_uri(taxon_id: int | None) -> str | None:
    """identifiers.org taxonomy URI derived from an NCBI taxon id, or None."""
    if taxon_id is None:
        return None
    return f"https://identifiers.org/taxonomy:{taxon_id}"

router = APIRouter()


@router.get("/genomes", response_model=PaginatedResponse[GenomeResponse])
def list_genomes(
    digest: str | None = Query(None, description="The genome digest"),
    alias: str | None = Query(None, description="The genome alias"),
    search: SearchParams = Depends(get_search_params),
    pagination: PaginationParams = Depends(get_pagination_params),
    session: Session = Depends(get_db_session),
    rgc: Refgenie = Depends(get_refgenie),
):
    """
    List genomes with their aliases and asset counts.
    Supports both legacy query parameters (digest, alias) and new search functionality.

    Legacy parameters:
    - digest: Filter by exact genome digest
    - alias: Filter by exact genome alias

    Searchable columns: digest, description, species_name, common_name,
    assembly_source, assembly_accession, aliases

    Aliases are sourced from the mode-selected alias manager (``rgc.alias``) rather
    than the SQL ``alias`` table, so this works in local (store-backed) mode where
    the SQL alias table is empty as well as server (SQL-backed) mode.
    """
    # Query genomes directly; aliases are attached from the manager afterward.
    # No join/distinct: genomes with zero aliases must still appear.
    base_genome_query = select(Genome)

    # Apply legacy filters first
    if digest:
        base_genome_query = base_genome_query.where(Genome.digest == digest)
    if alias:
        try:
            digest_from_alias = rgc.alias.resolve(alias)
        except MissingAliasError:
            # No genome matches this alias -> empty page, not 404.
            return PaginatedResponse(
                items=[],
                pagination=PaginationMeta(
                    offset=pagination.offset,
                    limit=pagination.limit,
                    total=0,
                ),
            )
        base_genome_query = base_genome_query.where(Genome.digest == digest_from_alias)

    # Apply search functionality. Alias search cannot use a SQL join (aliases may
    # live in the store), so build one combined OR condition: scalar columns via
    # SQL, plus a materialized Genome.digest IN (...) for alias matches.
    searchable_columns = [
        "digest",
        "description",
        "species_name",
        "common_name",
        "assembly_source",
        "assembly_accession",
        "aliases",
    ]
    search_term = search.q.strip() if search.q else ""
    if search_term:
        _validate_search_fields(search.search_fields, searchable_columns)
        fields_to_search = search.search_fields or searchable_columns
        valid_fields = [f for f in fields_to_search if f in searchable_columns]

        conditions = []
        for col in (
            "digest",
            "description",
            "species_name",
            "common_name",
            "assembly_source",
            "assembly_accession",
        ):
            if col in valid_fields:
                conditions.append(
                    _build_search_condition(getattr(Genome, col), search_term, search.operator)
                )
        if "aliases" in valid_fields:
            matching_digests = alias_digests_matching(rgc.alias, search_term, search.operator)
            conditions.append(Genome.digest.in_(matching_digests))

        if conditions:
            base_genome_query = base_genome_query.where(or_(*conditions))

    # Use paginate_query for database-level pagination
    paginated_result = paginate_query(session, base_genome_query, pagination)

    # Convert each genome to GenomeResponse with all aliases from the manager
    genome_responses = []
    for genome in paginated_result.items:
        aliases = rgc.alias.get_for_genome(genome.digest)

        genome_response = GenomeResponse(
            digest=genome.digest,
            description=genome.description,
            aliases=aliases,
            asset_count=len(genome.asset_groups),
            species_name=genome.species_name,
            common_name=genome.common_name,
            taxon_id=genome.taxon_id,
            assembly_source=genome.assembly_source,
            assembly_accession=genome.assembly_accession,
        )
        genome_responses.append(genome_response)

    return PaginatedResponse(
        items=genome_responses,
        pagination=paginated_result.pagination,
    )


@router.get("/genomes/{digest}", response_model=GenomeDetailResponse)
def get_genome(
    digest: str,
    session: Session = Depends(get_db_session),
    rgc: Refgenie = Depends(get_refgenie),
):
    """Return a genome's core columns plus its FHR sidecar blob and taxon URI.

    The FHR blob is read from the RefgetStore exactly as ``get_alias`` does, so
    the UI detail page gets the whole metadata long tail in one request.
    """
    genome = session.exec(select(Genome).where(Genome.digest == digest)).unique().one_or_none()
    if not genome:
        raise HTTPException(status_code=404, detail=f"Genome {digest} not found")

    detail = GenomeDetailResponse.model_validate(genome, from_attributes=True)
    detail.taxon_uri = _taxon_uri(genome.taxon_id)
    fhr = rgc.store_router.get_fhr_metadata(digest, store_name=genome.store_name)
    detail.fhr = fhr.to_dict() if fhr is not None else None
    return detail


@router.get("/asset_groups", response_model=PaginatedResponse[AssetGroupPublic])
def list_asset_groups(
    genome_digest: str | None = Query(None, description="The genome digest"),
    asset_class: str | None = Query(None, description="The asset class"),
    asset_group_name: str | None = Query(None, description="The asset group name"),
    asset_group_id: int | None = Query(None, description="The asset group ID"),
    search: SearchParams = Depends(get_search_params),
    pagination: PaginationParams = Depends(get_pagination_params),
    session: Session = Depends(get_db_session),
):
    """
    List asset groups.
    Supports both legacy query parameters and new search functionality.

    Legacy parameters:
    - genome_digest: Filter by exact genome digest
    - asset_class: Filter by exact asset class name
    - asset_group_name: Filter by exact asset group name
    - asset_group_id: Filter by exact asset group ID

    Searchable columns: name
    """
    query = select(AssetGroup)

    # Apply legacy filters first
    if genome_digest:
        query = query.join(Genome).where(Genome.digest == genome_digest)
    if asset_class:
        query = query.join(AssetClass).where(AssetClass.name == asset_class)
    if asset_group_name:
        query = query.where(AssetGroup.name == asset_group_name)
    if asset_group_id:
        query = query.where(AssetGroup.id == asset_group_id)

    # Apply search functionality
    searchable_columns = ["name"]
    query = apply_search_to_query(query, search, searchable_columns, AssetGroup)

    return paginate_query(session, query, pagination)


@router.get("/asset_groups/{id}", response_model=AssetGroupPublic)
def get_asset_group(id: int, session: Session = Depends(get_db_session)):
    asset_group = session.exec(select(AssetGroup).where(AssetGroup.id == id)).unique().one_or_none()
    if not asset_group:
        raise HTTPException(status_code=404, detail=f"Asset group {id} not found")
    return asset_group


@router.get("/assets/{digest}", response_model=AssetResponse)
def get_asset(digest: str, session: Session = Depends(get_db_session)):
    # Eager-load asset_group -> asset_class so the resolved serving_modes property
    # serializes without triggering lazy loads. seek_keys feeds
    # AssetResponse.seek_keys, which the web UI's asset page lists.
    asset = (
        session.exec(
            select(Asset)
            .options(
                selectinload(Asset.asset_group).selectinload(AssetGroup.asset_class),
                selectinload(Asset.asset_names),
                selectinload(Asset.seek_keys),
            )
            .where(Asset.digest == digest)
        )
        .unique()
        .one_or_none()
    )
    if not asset:
        raise HTTPException(status_code=404, detail=f"Asset {digest} not found")
    return asset


@router.get("/assets", response_model=PaginatedResponse[AssetResponse])
def list_assets(
    name: str | None = Query(None, description="The tag name"),
    asset_group_name: str | None = Query(None, description="The asset group name"),
    genome_digest: str | None = Query(None, description="The genome id"),
    recipe_name: str | None = Query(
        None, description="The recipe name that was used to create the asset"
    ),
    asset_group_id: int | None = Query(
        None, description="The asset group ID to filter assets by"
    ),
    search: SearchParams = Depends(get_search_params),
    pagination: PaginationParams = Depends(get_pagination_params),
    session: Session = Depends(get_db_session),
):
    """
    List assets.
    Supports both legacy query parameters and new search functionality.

    Legacy parameters:
    - name: Filter by exact asset name
    - asset_group_name: Filter by exact asset group name
    - genome_digest: Filter by exact genome digest
    - recipe_name: Filter by exact recipe name
    - asset_group_id: Filter by exact asset group ID

    Searchable columns: name, digest, path
    """
    # Eager-load asset_group -> asset_class so the resolved serving_modes property
    # serializes without triggering lazy loads during response serialization.
    # asset_names feeds AssetResponse.names, seek_keys feeds
    # AssetResponse.seek_keys. Both are eager-loaded here because pydantic
    # serializes every field of the response model: without the eager load the
    # listing degrades into one extra query per asset, not into no query.
    query = select(Asset).options(
        selectinload(Asset.asset_group).selectinload(AssetGroup.asset_class),
        selectinload(Asset.asset_names),
        selectinload(Asset.seek_keys),
    )

    # Apply legacy filters first. Names live in the assetname table, so a
    # ``name`` filter resolves through it -- this is what lets a client pull by
    # a non-canonical name.
    if name:
        query = query.join(AssetName, AssetName.asset_digest == Asset.digest).where(
            AssetName.name == name
        )
    if asset_group_id:
        query = query.join(AssetGroup).where(AssetGroup.id == asset_group_id)
    if asset_group_name:
        query = query.join(AssetGroup).where(AssetGroup.name == asset_group_name)
    if genome_digest:
        query = query.join(AssetGroup).join(Genome).where(Genome.digest == genome_digest)
    if recipe_name:
        query = query.join(Recipe).where(Recipe.name == recipe_name)

    # Apply search functionality
    searchable_columns = ["name", "digest", "path"]
    query = apply_search_to_query(query, search, searchable_columns, Asset)

    return paginate_query(session, query, pagination)


@router.get("/assets/{asset_digest}/files")
def list_asset_files(asset_digest: str, session: Session = Depends(get_db_session)):
    """
    List files available for individual download from an asset.

    Returns ``{"asset_digest": ..., "files": [...]}`` built from the
    directory_contents of the StagedAsset record where mode="file". The
    ``/assets/{asset_digest}/files/{file_path}`` download sibling lives in
    `refgenie.server.routers.version4` (the dash app has no file downloads).
    """
    sa = helpers.get_staged_asset(session, asset_digest, "file")

    if not sa:
        raise HTTPException(
            status_code=404, detail=f"Asset {asset_digest} is not staged for file-level serving"
        )

    return {"asset_digest": asset_digest, "files": sa.directory_contents}


@router.get("/asset_classes", response_model=PaginatedResponse[AssetClassPublic])
def list_asset_classes(
    name: str | None = Query(None, description="The asset class name"),
    version: str | None = Query(None, description="The asset class version"),
    search: SearchParams = Depends(get_search_params),
    pagination: PaginationParams = Depends(get_pagination_params),
    session: Session = Depends(get_db_session),
):
    """
    List asset classes.
    Supports both legacy query parameters (name, version) and new search functionality.

    Legacy parameters:
    - name: Filter by exact asset class name
    - version: Filter by exact asset class version

    Searchable columns: name, version, description
    """
    query = select(AssetClass)

    # Apply legacy filters first
    if name:
        query = query.where(AssetClass.name == name)
    if version:
        query = query.where(AssetClass.version == version)

    # Apply search functionality
    searchable_columns = ["name", "version", "description"]
    query = apply_search_to_query(query, search, searchable_columns, AssetClass)

    return paginate_query(session, query, pagination)


@router.get("/asset_classes/{id}", response_model=AssetClassPublic)
def get_asset_class(id: int, session: Session = Depends(get_db_session)):
    asset_class = session.exec(select(AssetClass).where(AssetClass.id == id)).unique().one_or_none()
    if not asset_class:
        raise HTTPException(status_code=404, detail=f"Asset class {id} not found")
    return asset_class


@router.get("/recipes", response_model=PaginatedResponse[RecipePublic])
def list_recipes(
    output_asset_class: str | None = Query(None, description="The asset class"),
    name: str | None = Query(None, description="The recipe name"),
    version: str | None = Query(None, description="The recipe version"),
    search: SearchParams = Depends(get_search_params),
    pagination: PaginationParams = Depends(get_pagination_params),
    session: Session = Depends(get_db_session),
):
    """
    List recipes.
    Supports both legacy query parameters and new search functionality.

    Legacy parameters:
    - output_asset_class: Filter by exact output asset class name
    - name: Filter by exact recipe name
    - version: Filter by exact recipe version

    Searchable columns: name, version, description
    """
    query = select(Recipe)

    # Apply legacy filters first
    if output_asset_class:
        query = query.join(AssetClass).where(AssetClass.name == output_asset_class)
    if name:
        query = query.where(Recipe.name == name)
    if version:
        query = query.where(Recipe.version == version)

    # Apply search functionality
    searchable_columns = ["name", "version", "description"]
    query = apply_search_to_query(query, search, searchable_columns, Recipe)

    return paginate_query(session, query, pagination)


@router.get("/recipes/{id}", response_model=RecipePublic)
def get_recipe(id: int, session: Session = Depends(get_db_session)):
    query = select(Recipe).where(Recipe.id == id)
    recipe = session.exec(query).unique().one_or_none()
    if not recipe:
        raise HTTPException(status_code=404, detail=f"Recipe {id} not found")
    return recipe


@router.get("/configurations", response_model=PaginatedResponse[ConfigurationPublic])
def list_configurations(
    pagination: PaginationParams = Depends(get_pagination_params),
    session: Session = Depends(get_db_session),
):
    return paginate_query(session, select(Configuration), pagination)


@router.get("/configurations/{id}", response_model=ConfigurationPublic)
def get_configuration(id: int, session: Session = Depends(get_db_session)):
    configuration = (
        session.exec(select(Configuration).where(Configuration.id == id)).unique().one_or_none()
    )
    if not configuration:
        raise HTTPException(status_code=404, detail=f"Configuration {id} not found")
    return configuration


@router.get("/staged_assets", response_model=PaginatedResponse[StagedAssetPublic])
def list_staged_assets(
    asset_digest: str | None = Query(None, description="The asset digest"),
    mode: str | None = Query(None, description="The staging mode (file or archive)"),
    search: SearchParams = Depends(get_search_params),
    pagination: PaginationParams = Depends(get_pagination_params),
    session: Session = Depends(get_db_session),
):
    """
    List staged assets.
    Supports query parameters (asset_digest, mode) and search functionality.

    Searchable columns: asset_digest, mode
    """
    query = select(StagedAsset)

    if asset_digest:
        query = query.where(StagedAsset.asset_digest == asset_digest)
    if mode:
        query = query.where(StagedAsset.mode == mode)

    searchable_columns = ["asset_digest", "mode"]
    query = apply_search_to_query(query, search, searchable_columns, StagedAsset)

    return paginate_query(session, query, pagination)


@router.get("/staged_assets/{id}", response_model=StagedAssetPublic)
def get_staged_asset(id: int, session: Session = Depends(get_db_session)):
    staged = session.get(StagedAsset, id)
    if not staged:
        raise HTTPException(status_code=404, detail=f"Staged asset {id} not found")
    return staged


# New relationship endpoints for asset parents and children
@router.get("/relationships/{asset_digest}")
def get_asset_relationships(
    asset_digest: str,
    expand: bool = Query(
        False,
        description="If true, return the full asset objects instead of bare digests",
    ),
    session: Session = Depends(get_db_session),
):
    """
    Get parent and child relationships for a specific asset.

    Returns a dictionary with 'parents' and 'children' keys containing
    either asset digests (expand=false) or full asset objects (expand=true).
    """
    asset = (
        session.exec(
            select(Asset)
            .options(
                selectinload(Asset.parents).selectinload(Asset.asset_group),
                selectinload(Asset.children).selectinload(Asset.asset_group),
                selectinload(Asset.parents).selectinload(Asset.asset_names),
                selectinload(Asset.children).selectinload(Asset.asset_names),
            )
            .where(Asset.digest == asset_digest)
        )
        .unique()
        .one_or_none()
    )

    if not asset:
        raise HTTPException(status_code=404, detail=f"Asset {asset_digest} not found")

    if expand:
        # AssetResponse (not AssetPublic) is the wire shape for an asset: it
        # carries the resolved serving_modes, the asset class name and the
        # asset's other names, which a client needs to decide how to pull it.
        parents = [AssetResponse.model_validate(parent) for parent in asset.parents]
        children = [AssetResponse.model_validate(child) for child in asset.children]
    else:
        # For non-expanded, just return digests
        parents = (
            [parent.digest for parent in asset.parents if parent.digest] if asset.parents else []
        )
        children = (
            [child.digest for child in asset.children if child.digest] if asset.children else []
        )

    return {"asset_digest": asset_digest, "parents": parents, "children": children}


@router.get("/aliases", response_model=PaginatedResponse[AliasPublic])
def list_aliases(
    name: str | None = Query(None, description="The alias name"),
    genome_digest: str | None = Query(None, description="The genome digest"),
    search: SearchParams = Depends(get_search_params),
    pagination: PaginationParams = Depends(get_pagination_params),
    rgc: Refgenie = Depends(get_refgenie),
):
    """
    List aliases.
    Supports both legacy query parameters (name, genome_digest) and new search functionality.

    Legacy parameters:
    - name: Filter by exact alias name
    - genome_digest: Filter by exact genome digest

    Searchable columns: name

    Aliases are sourced from the mode-selected alias manager (``rgc.alias``) rather
    than the SQL ``alias`` table, so this works in local (store-backed) mode as well
    as server (SQL-backed) mode.
    """
    records = list(rgc.alias.list_all(genome_digest=genome_digest))

    # Apply legacy exact-name filter
    if name:
        records = [r for r in records if r.name == name]

    # Apply search functionality (name only)
    searchable_columns = ["name"]
    search_term = search.q.strip() if search.q else ""
    if search_term:
        _validate_search_fields(search.search_fields, searchable_columns)
        fields_to_search = search.search_fields or searchable_columns
        if "name" in fields_to_search:
            records = [
                r for r in records if name_matches_operator(r.name, search_term, search.operator)
            ]

    aliases = [AliasPublic(name=r.name, genome_digest=r.genome_digest) for r in records]

    return paginate_list(aliases, pagination)


@router.get("/aliases/{name}", response_model=AliasResponse)
def get_alias(name: str, rgc: Refgenie = Depends(get_refgenie)):
    """Resolve a genome alias to its digest, sequence collection, and provenance.

    Goes through the mode-selected alias manager rather than querying the SQL
    ``alias`` table, which is empty in local (store-backed) mode.
    """
    try:
        digest = rgc.alias.resolve(name)
    except MissingAliasError:
        raise HTTPException(status_code=404, detail=f"Alias '{name}' not found")

    # Get seqcol level 2, routing to the store that owns this genome.
    router = rgc.store_router
    collection = router.get_collection_level2(digest)
    if collection is None:
        raise HTTPException(status_code=404, detail=f"Collection data not found for {digest}")

    # Get FHR metadata (may be None). get_fhr_metadata returns a gtars
    # FhrMetadata object; the response field is a plain dict, so convert.
    _fhr = router.get_fhr_metadata(digest)
    fhr = _fhr.to_dict() if _fhr is not None else None

    return AliasResponse(
        alias=name,
        digest=digest,
        source="server",
        collection=collection,
        fhr=fhr,
    )
