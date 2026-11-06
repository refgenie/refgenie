"""The local actions API: the web equivalent of the CLI dispatch.

Pull, build and genome init are submitted as jobs (202 + ``JobRef``); the
other actions run in the request and return ``ActionResult``. Handlers are
``def`` so FastAPI runs the blocking manager calls on its threadpool.
"""

from collections.abc import Callable

from fastapi import APIRouter, Depends, HTTPException, Request
from fastapi.routing import APIRoute

from refgenie.core import Refgenie
from refgenie.db.tables import DataChannelType
from refgenie.exceptions import MissingGenomeError, RefgenieError
from refgenie.models import DataChannelSyncReport, GenomeAlias, GenomeDigest
from refgenie.server.dependencies import get_refgenie
from refgenie.server.errors import ErrorCode, ErrorResponse, WebError, action_error
from refgenie.server.jobs import (
    BuildJobParams,
    GenomeInitJobParams,
    JobManager,
    JobRef,
    PullJobParams,
    get_job_manager,
)
from refgenie.server.local.schemas import (
    ActionResult,
    AliasSetRequest,
    DataChannelAddRequest,
    PreflightResult,
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


class _ClassifiedRoute(APIRoute):
    """Answers a library error with the envelope ``classify_exception`` picks for it."""

    def get_route_handler(self) -> Callable:
        handler = super().get_route_handler()

        async def route(request: Request):
            try:
                return await handler(request)
            except (RefgenieError, ValueError, WebError) as exc:
                raise action_error(exc) from exc

        return route


#: Documents the envelope for every non-2xx in the OpenAPI schema.
_ERROR_RESPONSES = {"4XX": {"model": ErrorResponse}, "5XX": {"model": ErrorResponse}}

# On the router, so every route gets both: the anti-CSRF header and the bridge policy.
router = APIRouter(
    prefix="/actions",
    tags=["Actions"],
    dependencies=[Depends(require_action_header), Depends(require_action_origin)],
    responses=_ERROR_RESPONSES,
    route_class=_ClassifiedRoute,
)


# ---------------------------------------------------------------------------
# Long-running: submitted to the JobManager, 202 + JobRef
# ---------------------------------------------------------------------------


@router.post("/pull", status_code=202, response_model=JobRef)
def pull(
    req: PullJobParams,
    request: Request,
    jobs: JobManager = Depends(get_job_manager),
    rgc: Refgenie = Depends(get_refgenie),
) -> JobRef:
    """Queue a pull. The runner never prompts (``force_large=True``).

    From a bridge origin, ``server_url`` must be a subscribed server (403
    ``server_not_subscribed`` otherwise), since a bridge page cannot
    subscribe. Same-origin callers may pass any server.
    """
    if req.server_url and is_bridge_request(request):
        server_url = rgc.servers.find_subscription(req.server_url)
        if server_url is None:
            raise WebError(
                "This refgenie is not subscribed to that server, so a remote page cannot "
                "pull from it. Open your local refgenie to confirm the pull there, or "
                "subscribe to the server first.",
                ErrorCode.SERVER_NOT_SUBSCRIBED,
                403,
            )
        req = req.model_copy(update={"server_url": server_url})
    return jobs.submit_pull(req)


@router.post("/build", status_code=202, response_model=JobRef)
def build(req: BuildJobParams, jobs: JobManager = Depends(get_job_manager)) -> JobRef:
    """Queue a build. Builds run one at a time."""
    return jobs.submit_build(req)


@router.post("/build/preflight", response_model=PreflightResult)
def build_preflight(req: BuildJobParams, rgc: Refgenie = Depends(get_refgenie)) -> PreflightResult:
    """Check a build without starting it. Always 200; ``ok`` says whether it would start."""
    return rgc.build.preflight(
        recipe_name=req.recipe,
        genome_alias=GenomeAlias(req.genome),
        asset_group_name=req.asset_group,
        asset_name=req.asset,
        recipe_version=req.recipe_version,
        params=req.params.to_build_params() if req.params else None,
        stage=req.stage,
    )


@router.post("/genomes", status_code=202, response_model=JobRef)
def genome_init(req: GenomeInitJobParams, jobs: JobManager = Depends(get_job_manager)) -> JobRef:
    """Queue a genome init: load a FASTA into the RefgetStore and, by default, build fasta."""
    return jobs.submit_genome_init(req)


# ---------------------------------------------------------------------------
# Synchronous curation: request/response, 200 + ActionResult
# ---------------------------------------------------------------------------


@router.delete("/assets/{asset_digest}", response_model=ActionResult)
def delete_asset(asset_digest: str, rgc: Refgenie = Depends(get_refgenie)) -> ActionResult:
    """Delete one asset by digest. 409 if other assets were built from it."""
    registry_path = rgc.asset.remove_by_digest(digest=asset_digest)
    return ActionResult(
        message=f"Asset {registry_path} deleted", data={"registry_path": registry_path}
    )


@router.delete("/genomes/{genome_digest}", response_model=ActionResult)
def delete_genome(genome_digest: str, rgc: Refgenie = Depends(get_refgenie)) -> ActionResult:
    """Delete a genome and all its aliases, asset groups and assets. The UI confirms."""
    digest = _genome_digest(genome_digest)
    rgc.genome.remove(digest)
    return ActionResult(message=f"Genome {digest} deleted", data={"genome_digest": digest})


@router.post("/aliases", response_model=ActionResult)
def set_alias(req: AliasSetRequest, rgc: Refgenie = Depends(get_refgenie)) -> ActionResult:
    """Register an alias for a genome that exists locally.

    The existence check is needed: ``set_genome_alias`` creates a placeholder
    genome for an unknown digest.
    """
    if not rgc.genome.exists(req.genome_digest):
        raise MissingGenomeError(genome=req.genome_digest)
    rgc.set_genome_alias(
        alias_name=GenomeAlias(req.alias), genome_digest=GenomeDigest(req.genome_digest)
    )
    return ActionResult(
        message=f"Alias '{req.alias}' set for genome {req.genome_digest}",
        data={"alias": req.alias, "genome_digest": req.genome_digest},
    )


@router.delete("/aliases/{alias_name}", response_model=ActionResult)
def delete_alias(alias_name: str, rgc: Refgenie = Depends(get_refgenie)) -> ActionResult:
    """Remove one alias (the genome stays)."""
    rgc.alias.remove(GenomeAlias(alias_name))
    return ActionResult(message=f"Alias '{alias_name}' removed", data={"alias": alias_name})


@router.post("/subscriptions", response_model=ActionResult)
def subscribe(req: SubscribeRequest, rgc: Refgenie = Depends(get_refgenie)) -> ActionResult:
    """Subscribe to servers; ``reset`` replaces the list. Returns the new list."""
    rgc.servers.subscribe(server_urls=req.server_urls, reset=req.reset)
    return ActionResult(
        message=f"Subscribed to: {', '.join(req.server_urls)}",
        data={"subscriptions": rgc.servers.subscriptions()},
    )


@router.delete("/subscriptions", response_model=ActionResult)
def unsubscribe(req: UnsubscribeRequest, rgc: Refgenie = Depends(get_refgenie)) -> ActionResult:
    """Unsubscribe from servers. Unknown URLs are ignored."""
    rgc.servers.unsubscribe(server_urls=req.server_urls)
    return ActionResult(
        message=f"Unsubscribed from: {', '.join(req.server_urls)}",
        data={"subscriptions": rgc.servers.subscriptions()},
    )


#: Appended to the add response. The CLI logs its own wording.
DATA_CHANNEL_WARNING = (
    "Data channels are not verified. Recipes from a channel run shell commands on "
    "this machine when you build; only add channels you trust."
)


def _sync_channel(rgc: Refgenie, name: str) -> DataChannelSyncReport:
    """Sync a channel, skipping items already present. 502 if it fails."""
    try:
        return rgc.sources.sync_channel(name, exists_ok=True)
    except RefgenieError as exc:  # channel missing or unreachable
        raise WebError(str(exc), ErrorCode.REFGENIE_ERROR, 502) from exc


@router.post("/data_channels", response_model=ActionResult)
def add_data_channel(
    req: DataChannelAddRequest, rgc: Refgenie = Depends(get_refgenie)
) -> ActionResult:
    """Record a data channel and, by default, sync it. The message carries the trust warning."""
    rgc.sources.add_channel(
        name=req.name,
        type=DataChannelType(req.channel_type),
        index_address=req.index_address,
        description=req.description,
    )
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
        raise HTTPException(status_code=404, detail=f"Data channel '{name}' not found.")
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
    """Make one asset its group's default."""
    rgc.asset.group.set_default(
        asset_group_name=req.asset_group,
        asset_name=req.asset,
        genome_digest=_genome_digest(req.genome_digest),
    )
    return ActionResult(
        message=f"Default asset for {req.genome_digest}/{req.asset_group} set to '{req.asset}'",
        data={
            "genome_digest": req.genome_digest,
            "asset_group": req.asset_group,
            "asset": req.asset,
        },
    )


def _genome_digest(value: str) -> GenomeDigest:
    """``value`` as a genome digest; a malformed one is a missing genome."""
    try:
        return GenomeDigest(value)
    except ValueError:
        raise MissingGenomeError(genome=value) from None
