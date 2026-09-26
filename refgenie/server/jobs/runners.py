"""What a job actually does: one function per kind, calling the root.

The work happens on a worker thread, and the answer arrives as a job record
instead of a response body. A failed job still says *which* kind of failure it
was (asset exists, missing genome, and so on), so the UI can tell the user what
to do next instead of showing a generic error.

The classification itself is not implemented here. `classify_exception` in
`refgenie/server/errors.py` is the single mapping, so an `AssetExistsError`
reports `asset_exists` whether it comes back from an HTTP call or off a job
record. Runners raise; the manager catches and classifies.
"""

from pathlib import Path

from refgenie.db.tables import Asset
from refgenie.logger import logger
from refgenie.models import BuildParams, GenomeAlias, GenomeDigest
from refgenie.server.errors import ErrorCode, WebError
from refgenie.server.jobs.manager import JobContext
from refgenie.server.jobs.schemas import (
    BuildJobParams,
    JobResult,
    BuildParamsSpec,
    GenomeInitJobParams,
    JobKind,
    PullJobParams,
)
from refgenie.server.jobs.tailer import build_log_tailer

__all__ = [
    "DEFAULT_RUNNERS",
    "build_runner",
    "genome_init_runner",
    "pull_runner",
    "to_build_params",
]


def to_build_params(spec: BuildParamsSpec | None) -> BuildParams | None:
    """The library-side `BuildParams` for a job's JSON-plain params spec."""
    if spec is None:
        return None
    return BuildParams(
        assets=spec.assets,
        params=spec.params,
        files={name: Path(path) for name, path in spec.files.items()} if spec.files else None,
    )


def _no_prompt(message: str) -> bool:
    """Refuse every confirmation, without reading stdin.

    Passing this explicitly matters. `resolve_confirmer(None)` falls through to
    an interactive `rich` prompt whenever `enable_interactive_prompts()` has
    been called in this process -- and it has been, because `refgenie dash` is
    launched from `main_cli`. Without an explicit confirmer a pull job would
    block a worker thread forever on a prompt against a terminal nobody is
    watching.
    """
    logger.info(f"Declined automatically (background job): {message}")
    return False


def pull_runner(ctx: JobContext) -> JobResult | None:
    """Pull one asset from a subscribed (or explicitly named) server."""
    params: PullJobParams = ctx.params
    # The wire keeps two fields; the library takes one typed identifier. A name
    # wins when both are set.
    genome = (
        GenomeAlias(params.genome_name)
        if params.genome_name
        else GenomeDigest(params.genome_digest)
    )
    target = f"{genome}/{params.asset_group_name}"
    # The phase names come from the vocabulary the UI renders
    # (frontend/src/components/jobs/phases.ts). Only the phases this layer can
    # observe truthfully are emitted here; the rest of a pull's phases
    # (`query`, `verify`, `unpack`, ...) are announced from inside the puller
    # and the download client, which are the only places that know them.
    ctx.phase("resolve", f"Pulling {target}")

    asset = ctx.refgenie.transfer.pull(
        asset_group_name=params.asset_group_name,
        genome=genome,
        asset_name=params.asset_name,
        force=params.force,
        # A browser cannot answer a "this is 12 GB, continue?" prompt, and the
        # user already clicked the button that means yes.
        force_large=True,
        force_server_urls=[params.server_url] if params.server_url else None,
        confirm=_no_prompt,
    )
    if asset is None:
        # transfer.pull() returns None only when there is nothing subscribed and the
        # subscribe prompt was declined -- which, with `_no_prompt`, is always.
        raise WebError(
            "No subscribed servers; subscribe to a server first.",
            ErrorCode.NO_SUBSCRIPTIONS,
            409,
        )
    return _asset_result(ctx, asset)


def build_runner(ctx: JobContext) -> JobResult | None:
    """Build one asset from a recipe, streaming pypiper's log to the browser."""
    params: BuildJobParams = ctx.params
    ctx.phase("resolve_recipe", f"Building {params.genome_name}/{params.asset_group_name}")

    # Build where the recipe says its tools live. Preflight already read the
    # asset name inside the recipe's image (`build.resolve_custom_seek_keys`), so a
    # container-based recipe that preflighted clean must not then run on a
    # host that has no bowtie2. The CLI keeps its explicit `--docker` flag;
    # the web has no flag, so the recipe decides.
    recipe = ctx.refgenie.recipe.get(params.recipe_name, params.recipe_version)
    docker = recipe.docker_image is not None

    with build_log_tailer(ctx, params.genome_name, params.asset_group_name):
        asset = ctx.refgenie.build.run(
            recipe_name=params.recipe_name,
            genome_alias=GenomeAlias(params.genome_name),
            asset_group_name=params.asset_group_name,
            asset_name=params.asset_name,
            recipe_version=params.recipe_version,
            asset_description=params.asset_description,
            stage=params.stage,
            pull_parents=params.pull_parents,
            params=to_build_params(params.params),
            docker=docker,
        )
    if asset is None:
        raise WebError(
            "Build failed. See the pipeline log for details.",
            ErrorCode.BUILD_FAILED,
            500,
        )
    return _asset_result(ctx, asset, staged=params.stage)


def genome_init_runner(ctx: JobContext) -> JobResult | None:
    """Initialize a genome from a FASTA, and build its fasta asset.

    Runs on the build executor: it ends in a build, and a build must never
    overlap another one.
    """
    params: GenomeInitJobParams = ctx.params
    ctx.phase("resolve", f"Initializing {params.aliases[0]}")

    with build_log_tailer(ctx, params.aliases[0], "fasta"):
        digest, created = ctx.refgenie.build.initialize_and_build(
            fasta_file_path=Path(params.fasta),
            genome_names=[GenomeAlias(a) for a in params.aliases],
            description=params.description or "",
            species_name=params.species,
            build_fasta=params.build_fasta_asset,
        )
    return JobResult(genome_digest=digest, created=created, aliases=params.aliases)


#: The runner for each kind. The manager takes a copy at construction, so a
#: test can substitute fakes without touching this table.
DEFAULT_RUNNERS = {
    JobKind.PULL: pull_runner,
    JobKind.BUILD: build_runner,
    JobKind.GENOME_INIT: genome_init_runner,
}


def _asset_result(ctx: JobContext, asset: Asset, staged: bool = False) -> JobResult:
    """The result payload for a pull or build.

    Enough for the frontend to refresh exactly the rows that changed instead of
    reloading the page. Read defensively: the returned `Asset` may be detached
    from the session that produced it, so any relationship access (which is
    what `registry_path` is) can raise -- and a successful pull must not be
    reported as a failure because its receipt could not be formatted.
    """
    result = JobResult(
        asset_digest=getattr(asset, "digest", None),
        asset_name=getattr(asset, "name", None),
        staged=staged,
    )
    try:
        result.registry_path = asset.registry_path
        result.asset_group_name = asset.asset_group.name
        result.genome_digest = asset.asset_group.genome.digest
    except Exception:  # noqa: BLE001 - a detached instance, not a failed job
        logger.debug(f"Job {ctx.job_id}: could not resolve the asset's registry path")
    return result
