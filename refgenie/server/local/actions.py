"""The local actions API: the web equivalent of the CLI dispatch.

Thin handlers over the ``Refgenie`` facade and the JobManager. Long-running
operations (pull, build, genome init) are **submitted**, never run in the
request: the handler translates its web request model into typed job params
and returns 202 with the ``JobRef``. Synchronous curation (aliases, deletes,
subscriptions, defaults) stays request/response and returns ``ActionResult``.

Every route requires the ``X-Refgenie-Action`` header (attached once, on the
router, so a newly added route cannot forget it) and every failure uses the
envelope from :mod:`refgenie.server.errors`. There are **no** GET routes here:
reads are ``/v4``, job polling is ``/v1/jobs``.

All handlers are ``def``, not ``async def``, so FastAPI runs the blocking
manager calls on its threadpool.
"""

from pathlib import Path

from fastapi import APIRouter, Depends, HTTPException

from refgenie.core import Refgenie
from refgenie.exceptions import (
    MissingAliasError,
    MissingAssetError,
    MissingAssetGroupError,
    MissingGenomeError,
)
from refgenie.models import BuildParams
from refgenie.server.dependencies import get_refgenie
from refgenie.server.errors import ErrorCode, ErrorResponse, action_error
from refgenie.server.jobs import (
    BuildJobParams,
    GenomeInitJobParams,
    JobManager,
    JobRef,
    PullJobParams,
    get_job_manager,
)
from refgenie.server.jobs.schemas import BuildParamsSpec
from refgenie.server.local.schemas import (
    ActionResult,
    AliasSetRequest,
    BuildRequest,
    GenomeInitRequest,
    PreflightResult,
    PullRequest,
    SetDefaultAssetRequest,
    SubscribeRequest,
    UnsubscribeRequest,
)
from refgenie.server.local.security import require_action_header, require_action_origin

__all__ = ["router"]

#: Documents the envelope for every non-2xx in the OpenAPI schema.
_ERROR_RESPONSES = {"4XX": {"model": ErrorResponse}, "5XX": {"model": ErrorResponse}}

# require_action_header is the anti-CSRF control (forces a preflight);
# require_action_origin is the bridge policy on top (which cross-origin
# callers may reach which action). Both attach once, on the router, so a
# newly added route cannot forget either.
router = APIRouter(
    prefix="/actions",
    tags=["Actions"],
    dependencies=[Depends(require_action_header), Depends(require_action_origin)],
    responses=_ERROR_RESPONSES,
)


# ---------------------------------------------------------------------------
# Long-running: submitted to the JobManager, 202 + JobRef
# ---------------------------------------------------------------------------


@router.post("/pull", status_code=202, response_model=JobRef)
def pull(req: PullRequest, jobs: JobManager = Depends(get_job_manager)) -> JobRef:
    """Queue a pull. The runner supplies ``force_large=True`` and a refusing
    confirmer -- a browser cannot answer prompts."""
    return jobs.submit_pull(
        PullJobParams(
            asset_group_name=req.asset_group,
            genome_name=req.genome,
            genome_digest=req.genome_digest,
            asset_name=req.asset,
            server_url=req.server_url,
            force=req.force,
        )
    )


@router.post("/build", status_code=202, response_model=JobRef)
def build(req: BuildRequest, jobs: JobManager = Depends(get_job_manager)) -> JobRef:
    """Queue a build. One at a time -- the build executor has one slot."""
    return jobs.submit_build(_build_job_params(req))


@router.post("/build/preflight", response_model=PreflightResult)
def build_preflight(req: BuildRequest, rgc: Refgenie = Depends(get_refgenie)) -> PreflightResult:
    """Validate a prospective build without starting it. Always 200; ``ok``
    says whether the build would start, ``errors`` is field-scoped."""
    return rgc.preflight_build(
        recipe_name=req.recipe,
        genome_name=req.genome,
        asset_group_name=req.asset_group,
        asset_name=req.asset,
        recipe_version=req.recipe_version,
        params=_to_build_params(req),
        stage=req.stage,
    )


@router.post("/genomes", status_code=202, response_model=JobRef)
def genome_init(req: GenomeInitRequest, jobs: JobManager = Depends(get_job_manager)) -> JobRef:
    """Queue a genome initialization: ingest a FASTA into the RefgetStore and
    (by default) build the fasta asset. Long-running, so it is a job."""
    return jobs.submit_genome_init(
        GenomeInitJobParams(
            fasta=req.fasta,
            aliases=req.aliases,
            description=req.description,
            species=req.species,
            build_fasta_asset=req.build_fasta_asset,
        )
    )


# ---------------------------------------------------------------------------
# Synchronous curation: request/response, 200 + ActionResult
# ---------------------------------------------------------------------------


@router.delete("/assets/{asset_digest}", response_model=ActionResult)
def delete_asset(asset_digest: str, rgc: Refgenie = Depends(get_refgenie)) -> ActionResult:
    """Delete one asset by digest. 409 (naming the children) if anything was
    built from it."""
    try:
        registry_path = rgc.asset.remove_by_digest(digest=asset_digest)
    except MissingAssetError:
        raise HTTPException(
            status_code=404,
            detail={
                "code": str(ErrorCode.ASSET_NOT_FOUND),
                "message": f"Asset identified with digest '{asset_digest}' not found.",
            },
        ) from None
    except ValueError as exc:  # has children; the message names them
        raise action_error(exc) from exc
    return ActionResult(
        message=f"Asset {registry_path} deleted", data={"registry_path": registry_path}
    )


@router.delete("/genomes/{genome_ref}", response_model=ActionResult)
def delete_genome(genome_ref: str, rgc: Refgenie = Depends(get_refgenie)) -> ActionResult:
    """Delete a genome (all aliases, asset groups and assets) by alias or
    digest. No ``force`` parameter: the CLI's force only suppresses a terminal
    prompt, and the web layer never prompts -- confirmation is the UI's job."""
    try:
        digest = rgc.alias.resolve(genome_ref)
    except MissingAliasError:
        digest = genome_ref
    try:
        rgc.genome.remove(digest)
    except MissingGenomeError as exc:
        raise action_error(exc) from exc
    return ActionResult(message=f"Genome {digest} deleted", data={"genome_digest": digest})


@router.post("/aliases", response_model=ActionResult)
def set_alias(req: AliasSetRequest, rgc: Refgenie = Depends(get_refgenie)) -> ActionResult:
    """Register an alias for a genome that exists locally.

    The existence guard matches the CLI and is not redundant:
    ``set_genome_alias`` *creates* a placeholder genome row for any unknown
    digest, so without the guard this endpoint would silently manufacture
    phantom genomes (the old dash endpoint's 404 branch could never fire).
    """
    if not rgc.genome.exists(req.genome_digest):
        raise action_error(MissingGenomeError(genome=req.genome_digest))
    rgc.set_genome_alias(alias_name=req.alias, genome_digest=req.genome_digest)
    return ActionResult(
        message=f"Alias '{req.alias}' set for genome {req.genome_digest}",
        data={"alias": req.alias, "genome_digest": req.genome_digest},
    )


@router.delete("/aliases/{alias_name}", response_model=ActionResult)
def delete_alias(alias_name: str, rgc: Refgenie = Depends(get_refgenie)) -> ActionResult:
    """Remove one alias (the genome stays)."""
    try:
        rgc.alias.remove(alias_name)
    except MissingAliasError as exc:
        raise action_error(exc) from exc
    return ActionResult(message=f"Alias '{alias_name}' removed", data={"alias": alias_name})


@router.post("/subscriptions", response_model=ActionResult)
def subscribe(req: SubscribeRequest, rgc: Refgenie = Depends(get_refgenie)) -> ActionResult:
    """Subscribe to one or more servers; ``reset`` replaces instead of adding.
    ``data.subscriptions`` is the resulting list, so the UI needs no second
    request to re-render."""
    rgc.configuration.subscribe(server_urls=req.server_urls, reset=req.reset)
    return ActionResult(
        message=f"Subscribed to: {', '.join(req.server_urls)}",
        data={"subscriptions": list(rgc.configuration.get_server_subscriptions())},
    )


@router.delete("/subscriptions", response_model=ActionResult)
def unsubscribe(req: UnsubscribeRequest, rgc: Refgenie = Depends(get_refgenie)) -> ActionResult:
    """Unsubscribe from one or more servers. Unknown URLs are a no-op, matching
    the manager's set semantics."""
    rgc.configuration.unsubscribe(server_urls=req.server_urls)
    return ActionResult(
        message=f"Unsubscribed from: {', '.join(req.server_urls)}",
        data={"subscriptions": list(rgc.configuration.get_server_subscriptions())},
    )


@router.post("/assets/default", response_model=ActionResult)
def set_default_asset(
    req: SetDefaultAssetRequest, rgc: Refgenie = Depends(get_refgenie)
) -> ActionResult:
    """Flag one asset as its group's default (single-transaction swap)."""
    try:
        rgc.asset.set_default(
            asset_group_name=req.asset_group,
            asset_name=req.asset,
            genome_digest=req.genome_digest,
        )
    except (MissingAssetError, MissingAssetGroupError, MissingGenomeError) as exc:
        raise action_error(exc) from exc
    return ActionResult(
        message=f"Default asset for {req.genome_digest}/{req.asset_group} set to '{req.asset}'",
        data={
            "genome_digest": req.genome_digest,
            "asset_group": req.asset_group,
            "asset": req.asset,
        },
    )


# ---------------------------------------------------------------------------
# Mapping helpers
# ---------------------------------------------------------------------------


def _build_job_params(req: BuildRequest) -> BuildJobParams:
    """Web field names -> facade-flavored job params."""
    return BuildJobParams(
        recipe_name=req.recipe,
        genome_name=req.genome,
        asset_group_name=req.asset_group,
        asset_name=req.asset,
        recipe_version=req.recipe_version,
        asset_description=req.description,
        stage=req.stage,
        pull_parents=req.pull_parents,
        params=BuildParamsSpec(**req.params.model_dump()) if req.params else None,
    )


def _to_build_params(req: BuildRequest) -> BuildParams | None:
    """The library-side ``BuildParams`` for a preflight call."""
    if req.params is None:
        return None
    return BuildParams(
        assets=req.params.assets,
        params=req.params.params,
        files=(
            {name: Path(path) for name, path in req.params.files.items()}
            if req.params.files
            else None
        ),
    )
