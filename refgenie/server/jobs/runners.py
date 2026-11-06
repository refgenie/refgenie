"""One runner function per job kind.

Runners raise on failure; the manager catches and classifies the exception.
"""

from pathlib import Path

from refgenie.db.tables import Asset
from refgenie.exceptions import PullSkipped
from refgenie.logger import logger
from refgenie.models import GenomeAlias, GenomeDigest
from refgenie.server.errors import ErrorCode, WebError
from refgenie.server.jobs.manager import JobContext
from refgenie.server.jobs.schemas import (
    BuildJobParams,
    GenomeInitJobParams,
    JobResult,
    JobKind,
    PullJobParams,
)
from refgenie.server.jobs.tailer import build_log_tailer

__all__ = [
    "DEFAULT_RUNNERS",
    "build_runner",
    "genome_init_runner",
    "pull_runner",
]


def _no_prompt(message: str) -> bool:
    """Refuse every confirmation, without reading stdin.

    Without it, a job under `refgenie dash` would block on an interactive prompt.
    """
    logger.info(f"Declined automatically (background job): {message}")
    return False


def pull_runner(ctx: JobContext) -> JobResult | None:
    """Pull one asset from a subscribed (or explicitly named) server."""
    params: PullJobParams = ctx.params
    genome = GenomeAlias(params.genome) if params.genome else GenomeDigest(params.genome_digest)
    target = f"{genome}/{params.asset_group}"
    # Later phases are announced from inside the puller and download client.
    ctx.phase("resolve", f"Pulling {target}")

    asset = ctx.refgenie.transfer.pull(
        asset_group_name=params.asset_group,
        genome=genome,
        asset_name=params.asset,
        force=params.force,
        # The user already confirmed by clicking Pull.
        force_large=True,
        force_server_urls=[params.server_url] if params.server_url else None,
        confirm=_no_prompt,
    )
    # Failures raise; None means the pull was skipped.
    if asset is None:
        raise PullSkipped(f"Skipped pull of {target}: existing content was kept.")
    return _asset_result(ctx, asset)


def build_runner(ctx: JobContext) -> JobResult | None:
    """Build one asset from a recipe, streaming pypiper's log to the browser."""
    params: BuildJobParams = ctx.params
    ctx.phase("resolve_recipe", f"Building {params.genome}/{params.asset_group}")

    # The web has no `--container` flag: a recipe with an image builds in it.
    recipe = ctx.refgenie.recipe.get(params.recipe, params.recipe_version)
    container = recipe.docker_image is not None

    with build_log_tailer(ctx, params.genome, params.asset_group):
        asset = ctx.refgenie.build.run(
            recipe_name=params.recipe,
            genome_alias=GenomeAlias(params.genome),
            asset_group_name=params.asset_group,
            asset_name=params.asset,
            recipe_version=params.recipe_version,
            asset_description=params.description,
            stage=params.stage,
            pull_parents=params.pull_parents,
            params=params.params.to_build_params() if params.params else None,
            container=container,
        )
    if asset is None:
        raise WebError(
            "Build failed. See the pipeline log for details.",
            ErrorCode.BUILD_FAILED,
            500,
        )
    return _asset_result(ctx, asset, staged=params.stage)


def genome_init_runner(ctx: JobContext) -> JobResult | None:
    """Initialize a genome from a FASTA, and build its fasta asset."""
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


#: The runner for each kind; the manager copies it at construction.
DEFAULT_RUNNERS = {
    JobKind.PULL: pull_runner,
    JobKind.BUILD: build_runner,
    JobKind.GENOME_INIT: genome_init_runner,
}


def _asset_result(ctx: JobContext, asset: Asset, staged: bool = False) -> JobResult:
    """The result payload for a pull or build.

    Relationship access may fail on a detached `Asset`; that is logged, not a
    failed job.
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
