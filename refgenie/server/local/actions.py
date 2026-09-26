"""The local actions API: the web equivalent of the CLI dispatch.

Thin handlers over the ``Refgenie`` root and the JobManager. Long-running
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

from fastapi import APIRouter, Depends, HTTPException, Request

from refgenie.core import Refgenie
from refgenie.db.tables import DataChannelType
from refgenie.exceptions import (
    MissingAliasError,
    MissingAssetError,
    MissingAssetGroupError,
    MissingGenomeError,
    RefgenieError,
)
from refgenie.models import BuildParams, DataChannelSyncReport, GenomeAlias, GenomeDigest
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
    DataChannelAddRequest,
    GenomeInitRequest,
    PreflightResult,
    PullRequest,
    SetDefaultAssetRequest,
    SubscribeRequest,
    UnsubscribeRequest,
)
from refgenie.server.local.security import (
    is_bridge_request,
    require_action_header,
    require_action_origin,
)

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
def pull(
    req: PullRequest,
    request: Request,
    jobs: JobManager = Depends(get_job_manager),
    rgc: Refgenie = Depends(get_refgenie),
) -> JobRef:
    """Queue a pull. The runner supplies ``force_large=True`` and a refusing
    confirmer -- a browser cannot answer prompts.

    A ``server_url`` from a bridge origin must name a subscription (403
    ``server_not_subscribed`` otherwise), because a bridge page cannot
    subscribe, so the subscription list is the only place the user said which
    servers they trust. A same-origin caller passes any ``server_url``: it
    could subscribe first anyway, and the local ``/pull`` confirmation page is
    exactly how a user pulls from a server they are not subscribed to.
    """
    server_url = req.server_url
    if server_url and is_bridge_request(request):
        server_url = rgc.servers.find_subscription(server_url)
        if server_url is None:
            raise HTTPException(
                status_code=403,
                detail={
                    "code": str(ErrorCode.SERVER_NOT_SUBSCRIBED),
                    "message": (
                        "This refgenie is not subscribed to that server, so a remote "
                        "page cannot pull from it. Open your local refgenie to confirm "
                        "the pull there, or subscribe to the server first."
                    ),
                },
            )
    return jobs.submit_pull(
        PullJobParams(
            asset_group_name=req.asset_group,
            genome_name=req.genome,
            genome_digest=req.genome_digest,
            asset_name=req.asset,
            server_url=server_url,
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
    return rgc.build.preflight(
        recipe_name=req.recipe,
        genome_alias=GenomeAlias(req.genome),
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


@router.delete("/genomes/{genome_digest}", response_model=ActionResult)
def delete_genome(genome_digest: str, rgc: Refgenie = Depends(get_refgenie)) -> ActionResult:
    """Delete a genome (all aliases, asset groups and assets) by digest. No
    ``force`` parameter: the CLI's force only suppresses a terminal prompt, and
    the web layer never prompts -- confirmation is the UI's job."""
    digest = _genome_digest(genome_digest)
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
    phantom genomes and its 404 branch could never fire.
    """
    if not rgc.genome.exists(req.genome_digest):
        raise action_error(MissingGenomeError(genome=req.genome_digest))
    try:
        alias = GenomeAlias(req.alias)
    except ValueError as exc:
        raise action_error(exc) from exc
    rgc.set_genome_alias(alias_name=alias, genome_digest=GenomeDigest(req.genome_digest))
    return ActionResult(
        message=f"Alias '{req.alias}' set for genome {req.genome_digest}",
        data={"alias": req.alias, "genome_digest": req.genome_digest},
    )


@router.delete("/aliases/{alias_name}", response_model=ActionResult)
def delete_alias(alias_name: str, rgc: Refgenie = Depends(get_refgenie)) -> ActionResult:
    """Remove one alias (the genome stays)."""
    try:
        rgc.alias.remove(GenomeAlias(alias_name))
    except MissingAliasError as exc:
        raise action_error(exc) from exc
    return ActionResult(message=f"Alias '{alias_name}' removed", data={"alias": alias_name})


@router.post("/subscriptions", response_model=ActionResult)
def subscribe(req: SubscribeRequest, rgc: Refgenie = Depends(get_refgenie)) -> ActionResult:
    """Subscribe to one or more servers; ``reset`` replaces instead of adding.
    ``data.subscriptions`` is the resulting list, so the UI needs no second
    request to re-render."""
    rgc.servers.subscribe(server_urls=req.server_urls, reset=req.reset)
    return ActionResult(
        message=f"Subscribed to: {', '.join(req.server_urls)}",
        data={"subscriptions": rgc.servers.subscriptions()},
    )


@router.delete("/subscriptions", response_model=ActionResult)
def unsubscribe(req: UnsubscribeRequest, rgc: Refgenie = Depends(get_refgenie)) -> ActionResult:
    """Unsubscribe from one or more servers. Unknown URLs are a no-op, matching
    the manager's set semantics."""
    rgc.servers.unsubscribe(server_urls=req.server_urls)
    return ActionResult(
        message=f"Unsubscribed from: {', '.join(req.server_urls)}",
        data={"subscriptions": rgc.servers.subscriptions()},
    )


#: Shown on add from both surfaces; the CLI logs the same warning.
DATA_CHANNEL_WARNING = (
    "Data channels are not verified. Recipes from a channel run shell commands on "
    "this machine when you build; only add channels you trust."
)


def _sync_channel(rgc: Refgenie, name: str) -> DataChannelSyncReport:
    """Sync with ``exists_ok``: re-syncing a channel from the UI is a refresh,
    not a request to fail on everything it already has."""
    try:
        return rgc.sources.sync_channel(name, exists_ok=True)
    except RefgenieError as exc:  # channel missing or unreachable
        raise HTTPException(
            status_code=502,
            detail={"code": str(ErrorCode.REFGENIE_ERROR), "message": str(exc)},
        ) from exc


@router.post("/data_channels", response_model=ActionResult)
def add_data_channel(
    req: DataChannelAddRequest, rgc: Refgenie = Depends(get_refgenie)
) -> ActionResult:
    """Record a data channel and, by default, sync it. Every channel is
    untrusted for now, so the message carries the warning rather than a
    separate acknowledgement step blocking the add."""
    try:
        rgc.sources.add_channel(
            name=req.name,
            type=DataChannelType(req.channel_type),
            index_address=req.index_address,
            description=req.description,
        )
    except ValueError as exc:  # duplicate name, or '__' in it
        raise action_error(exc) from exc
    data: dict = {"name": req.name, "trusted": False}
    if req.sync:
        data["sync"] = _sync_channel(rgc, req.name).model_dump()
    return ActionResult(
        message=f"Data channel '{req.name}' added. {DATA_CHANNEL_WARNING}", data=data
    )


@router.delete("/data_channels/{name}", response_model=ActionResult)
def remove_data_channel(name: str, rgc: Refgenie = Depends(get_refgenie)) -> ActionResult:
    """Forget a data channel. Recipes already synced from it stay registered."""
    if not rgc.sources.remove_channel(name):
        raise HTTPException(
            status_code=404,
            detail={
                "code": str(ErrorCode.NOT_FOUND),
                "message": f"Data channel '{name}' not found.",
            },
        )
    return ActionResult(message=f"Data channel '{name}' removed", data={"name": name})


@router.post("/data_channels/{name}/sync", response_model=ActionResult)
def sync_data_channel(name: str, rgc: Refgenie = Depends(get_refgenie)) -> ActionResult:
    """Fetch the channel's asset classes and recipes again, skipping known ones."""
    report = _sync_channel(rgc, name)
    added = report.asset_classes_added + report.recipes_added
    message = (
        f"Synced '{name}': {added} new item(s)"
        if report.ok
        else f"Synced '{name}' with {len(report.errors)} failure(s)"
    )
    return ActionResult(message=message, data={"sync": report.model_dump()})


@router.post("/assets/default", response_model=ActionResult)
def set_default_asset(
    req: SetDefaultAssetRequest, rgc: Refgenie = Depends(get_refgenie)
) -> ActionResult:
    """Flag one asset as its group's default (single-transaction swap)."""
    genome_digest = _genome_digest(req.genome_digest)
    try:
        rgc.asset.group.set_default(
            asset_group_name=req.asset_group,
            asset_name=req.asset,
            genome_digest=genome_digest,
        )
    except (
        MissingAssetError,
        MissingAssetGroupError,
        MissingGenomeError,
    ) as exc:
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


def _genome_digest(value: str) -> GenomeDigest:
    """``value`` as a genome digest. A malformed one is a genome that is not here."""
    try:
        return GenomeDigest(value)
    except ValueError:
        raise action_error(MissingGenomeError(genome=value)) from None


def _build_job_params(req: BuildRequest) -> BuildJobParams:
    """Web field names -> Refgenie-flavored job params."""
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
