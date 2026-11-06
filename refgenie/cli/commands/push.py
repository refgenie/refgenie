"""The `push` command: model and handler."""

from typing import Literal

from pydantic import AliasChoices, BaseModel, Field

from refgenie.cli.commands.helpers import genome_digest_from
from refgenie.cli.errors import EXIT_NOT_FOUND, fail
from refgenie.exceptions import MissingRemoteError
from refgenie.logger import logger


class PushModel(BaseModel):
    """push: upload staged assets to cloud remotes."""

    remote: str | None = Field(
        default=None,
        description="Push only to this remote (by name or id). Default: push to all remotes with unpushed assets.",
        validation_alias=AliasChoices("r", "remote"),
    )
    genome: str | None = Field(
        default=None,
        description="Push only assets for this genome alias.",
        validation_alias=AliasChoices("g", "genome"),
    )
    genome_digest: str | None = Field(
        default=None,
        description="Push only assets for this genome digest, in place of an alias.",
        validation_alias=AliasChoices("genome-digest"),
    )
    dry_run: bool = Field(
        default=False,
        description="Show what would be pushed without executing.",
        validation_alias=AliasChoices("n", "dry-run"),
    )
    strategy: Literal["per_asset", "folder_sync"] = Field(
        default="per_asset",
        description="Push strategy: 'per_asset' (upload each asset) or 'folder_sync' (sync entire genome_stage_folder).",
    )


def _finish_push(
    pushed_count: int, skipped_count: int, failed_count: int, total: int, dry_run: bool
) -> None:
    """Log the push summary (omitting zero categories) and apply the exit contract.

    Exits non-zero if any upload failed, or if nothing was pushed but
    something was skipped.
    """
    parts = []
    if pushed_count:
        parts.append(f"{pushed_count} pushed")
    if skipped_count:
        parts.append(f"{skipped_count} skipped (no push_command)")
    if failed_count:
        parts.append(f"{failed_count} failed")
    logger.info(f"Push complete: {', '.join(parts)}")
    if dry_run:
        logger.info("(dry run -- no uploads were executed)")
    if failed_count:
        fail(f"push: {failed_count} of {total} uploads failed")
    if pushed_count == 0 and skipped_count > 0:
        fail(f"push: nothing pushed ({skipped_count} skipped, no push_command)")


def handle_push(cmd, refgenie) -> None:
    """Push staged assets to cloud remotes; see ``refgenie.managers.push``."""
    from refgenie.managers.push import push

    genome_digest = None
    if cmd.genome or cmd.genome_digest:
        genome_digest = genome_digest_from(refgenie, cmd.genome, cmd.genome_digest)
    try:
        report = push(
            refgenie,
            remote=cmd.remote,
            genome_digest=genome_digest,
            strategy=cmd.strategy,
            dry_run=cmd.dry_run,
        )
    except MissingRemoteError as e:
        fail(f"push: {e}", EXIT_NOT_FOUND)
    if not report.total:
        logger.info("Nothing to push. All remote-asset links are up to date.")
        return
    _finish_push(report.pushed, report.skipped, report.failed, report.total, cmd.dry_run)
