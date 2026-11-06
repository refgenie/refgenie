"""Upload staged assets to their push remotes.

Each remote's ``push_command`` is a shell template. On an S3 prefix an upload
is verified with the ``aws`` CLI before its link is marked pushed, and an
archive already present is not uploaded again (archives are content-addressed).
"""

import json
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from refgenie.logger import logger
from refgenie.utils.staging import staged_archive_path

__all__ = ["PushReport", "push"]


@dataclass
class PushReport:
    """What a push did, per remote-asset link."""

    total: int = 0
    pushed: int = 0
    #: Links whose remote has no push_command.
    skipped: int = 0
    failed: int = 0


def push(
    rgc,
    remote: str | None = None,
    genome_digest: str | None = None,
    strategy: Literal["per_asset", "folder_sync"] = "per_asset",
    dry_run: bool = False,
) -> PushReport:
    """Push every unpushed staged asset, optionally for one remote or genome.

    ``per_asset`` runs the push_command once per asset; ``folder_sync`` runs it
    once per remote over the whole genome_stage_folder.

    Raises:
        MissingRemoteError: If ``remote`` names no remote.
        ValueError: If there is work to do but genome_stage_folder is not set.
    """
    results = rgc.remote.unpushed(remote, genome_digest=genome_digest)
    report = PushReport(total=len(results))
    if not results:
        return report
    stage_folder = rgc.genome_stage_folder
    if stage_folder is None:
        raise ValueError("genome_stage_folder is not configured")
    if strategy == "folder_sync":
        _push_folder_sync(rgc, results, Path(stage_folder), dry_run, report)
    else:
        _push_per_asset(rgc, results, Path(stage_folder), dry_run, report)
    return report


def _push_per_asset(rgc, results, stage_folder: Path, dry_run: bool, report: PushReport) -> None:
    for link, remote, staged_asset in results:
        tag = f"{link.asset_digest} mode={link.mode}"
        if remote.push_command is None:
            logger.warning(
                f"Remote '{remote.name}' (id={remote.id}) has no push_command. Skipping {tag}."
            )
            report.skipped += 1
            continue
        local_path = _local_path(staged_asset, stage_folder)
        if local_path is None:
            logger.error(f"Asset records for {link.asset_digest} not found in DB. Cannot push.")
            report.failed += 1
            continue
        relative_path = local_path.relative_to(stage_folder)
        s3 = _s3_key(remote.prefix, relative_path)

        if link.mode == "archive" and s3 is not None and _s3_head_object(*s3) is not None:
            url = f"s3://{s3[0]}/{s3[1]}"
            if dry_run:
                logger.info(f"[dry-run] Would skip (already present): {url}")
            else:
                rgc.remote.mark_pushed(remote.id, link.asset_digest, link.mode)
                logger.info(f"Skipped (already present): {url} [{link.asset_digest}]")
            report.pushed += 1
            continue

        command = _render(
            remote.push_command,
            local_path=local_path,
            relative_path=relative_path,
            prefix=remote.prefix,
        )
        if dry_run:
            logger.info(f"[dry-run] Would execute: {command}")
            logger.info(f"  asset_digest={link.asset_digest} mode={link.mode} -> {remote.name}")
            report.pushed += 1
            continue
        logger.info(f"Pushing {local_path} to {remote.name}: {command}")
        if not _run_push_cmd(command, f"Push failed for {link.asset_digest} to {remote.name}"):
            report.failed += 1
            continue
        if s3 is not None and not _verify_s3_object(
            *s3, staged_asset.tarball_size, is_dir=local_path.is_dir()
        ):
            logger.error(
                f"Not marking {tag} as pushed: object did not verify at "
                f"s3://{s3[0]}/{s3[1]} after upload."
            )
            report.failed += 1
            continue
        rgc.remote.mark_pushed(remote.id, link.asset_digest, link.mode)
        report.pushed += 1
        logger.info(f"Pushed {tag} to {remote.name}")


def _push_folder_sync(rgc, results, stage_folder: Path, dry_run: bool, report: PushReport) -> None:
    """One sync per remote; each link is then marked pushed only if its file verifies."""
    by_remote = defaultdict(list)
    for link, remote, staged_asset in results:
        by_remote[remote.id].append((link, remote, staged_asset))

    for items in by_remote.values():
        remote = items[0][1]
        if remote.push_command is None:
            logger.warning(f"Remote '{remote.name}' has no push_command. Skipping.")
            report.skipped += len(items)
            continue
        command = _render(
            remote.push_command, genome_stage_folder=stage_folder, prefix=remote.prefix
        )
        if dry_run:
            logger.info(f"[dry-run] Would execute: {command}")
            logger.info(f"  Would mark {len(items)} assets as pushed to {remote.name}")
            report.pushed += len(items)
            continue
        logger.info(f"Syncing genome_stage_folder to {remote.name}: {command}")
        if not _run_push_cmd(command, f"Folder sync failed for {remote.name}"):
            report.failed += len(items)
            continue

        marked = 0
        for link, _, staged_asset in items:
            tag = f"{link.asset_digest} mode={link.mode}"
            local_path = _local_path(staged_asset, stage_folder)
            if local_path is None or not local_path.exists():
                logger.error(
                    f"Not marking {tag} as pushed: staged file missing locally "
                    f"({local_path}), so the sync did not upload it."
                )
                report.failed += 1
                continue
            s3 = _s3_key(remote.prefix, local_path.relative_to(stage_folder))
            if s3 is not None and not _verify_s3_object(
                *s3, staged_asset.tarball_size, is_dir=local_path.is_dir()
            ):
                logger.error(
                    f"Not marking {tag} as pushed: object did not verify at "
                    f"s3://{s3[0]}/{s3[1]} after sync."
                )
                report.failed += 1
                continue
            rgc.remote.mark_pushed(remote.id, link.asset_digest, link.mode)
            marked += 1
        report.pushed += marked
        logger.info(f"Synced and marked {marked}/{len(items)} assets as pushed to {remote.name}")


def _local_path(staged_asset, stage_folder: Path) -> Path | None:
    """Where a staged asset sits on disk, or None if its asset rows are missing."""
    asset = staged_asset.asset
    asset_group = asset.asset_group if asset else None
    if not asset_group:
        return None
    if staged_asset.mode == "archive":
        return staged_archive_path(
            stage_folder, asset_group.genome_digest, asset_group.name, asset.digest
        )
    return stage_folder / asset_group.genome_digest / asset_group.name / asset.name


def _render(template: str, **values) -> str:
    """Fill ``{name}`` placeholders in a push_command template."""
    for name, value in values.items():
        template = template.replace(f"{{{name}}}", str(value))
    return template


def _run_push_cmd(command: str, failure: str) -> bool:
    """Run one push_command through the shell; log ``failure`` and stderr if it fails."""
    try:
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"{failure}: {e}")
        if e.stderr:
            logger.error(f"stderr: {e.stderr.strip()}")
        return False
    if result.stdout.strip():
        logger.debug(f"Push stdout: {result.stdout.strip()}")
    return True


def _parse_s3_prefix(prefix: str) -> tuple[str, str] | None:
    """``s3://bucket/path`` as ``(bucket, path)``; None if the prefix is not S3."""
    if not prefix.startswith("s3://"):
        return None
    bucket, _, path = prefix[5:].partition("/")
    return bucket, path.rstrip("/")


def _s3_key(prefix: str, relative_path: Path) -> tuple[str, str] | None:
    """``(bucket, key)`` for a path under an S3 prefix; None if the prefix is not S3."""
    parts = _parse_s3_prefix(prefix)
    if parts is None:
        return None
    bucket, base = parts
    return bucket, f"{base}/{relative_path}" if base else str(relative_path)


def _aws_json(*args: str) -> dict | None:
    """Run an ``aws s3api`` call and parse its JSON; None if it fails or aws is missing.

    ``--output json`` is explicit: a profile's ``output = text`` would make
    every object look absent.
    """
    try:
        result = subprocess.run(
            ["aws", "s3api", *args, "--output", "json"],
            check=True,
            capture_output=True,
            text=True,
        )
        return json.loads(result.stdout)
    except (subprocess.CalledProcessError, FileNotFoundError, json.JSONDecodeError):
        return None


def _s3_head_object(bucket: str, key: str) -> dict | None:
    """The object's metadata, or None if it is absent or the call fails."""
    return _aws_json("head-object", "--bucket", bucket, "--key", key)


def _s3_prefix_has_objects(bucket: str, prefix: str) -> bool:
    """Whether at least one object exists under ``s3://bucket/prefix/``."""
    listing = _aws_json(
        "list-objects-v2",
        "--bucket",
        bucket,
        "--prefix",
        f"{prefix.rstrip('/')}/",
        "--max-keys",
        "1",
    )
    return bool(listing and listing.get("Contents"))


def _verify_s3_object(
    bucket: str, key: str, expected_size: int | None, is_dir: bool = False
) -> bool:
    """Whether a push landed at ``s3://bucket/key``.

    A folder must have at least one object under ``key/``. A file's
    ``ContentLength`` must equal ``expected_size``; with none recorded
    (mode="file"), presence is enough.
    """
    if is_dir:
        if _s3_prefix_has_objects(bucket, key):
            return True
        logger.warning(f"Verification failed: s3://{bucket}/{key}/ is empty after push.")
        return False
    head = _s3_head_object(bucket, key)
    if head is None:
        logger.warning(f"Verification failed: s3://{bucket}/{key} not found after push.")
        return False
    if expected_size is None:
        logger.debug(f"No tarball_size recorded for s3://{bucket}/{key}; verifying presence only.")
        return True
    actual_size = head.get("ContentLength")
    if actual_size != expected_size:
        logger.warning(
            f"Verification failed: s3://{bucket}/{key} size {actual_size} "
            f"!= expected {expected_size}."
        )
        return False
    return True
