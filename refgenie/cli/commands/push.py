"""The `push` command: model and handler."""

import json
import subprocess
from pathlib import Path

from pydantic import AliasChoices, BaseModel, Field

from refgenie.cli.errors import fail
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
        description="Push only assets for this genome.",
        validation_alias=AliasChoices("g", "genome"),
    )
    dry_run: bool = Field(
        default=False,
        description="Show what would be pushed without executing.",
        validation_alias=AliasChoices("n", "dry-run"),
    )
    strategy: str = Field(
        default="per_asset",
        description="Push strategy: 'per_asset' (upload each asset) or 'folder_sync' (sync entire genome_stage_folder).",
    )


def _resolve_local_path(session, link, genome_stage_folder: Path) -> Path | None:
    """Resolve the staged local path for a RemoteAssetLink.

    Returns None if the link's Asset or AssetGroup rows are missing from the DB.
    """
    from sqlmodel import select

    from refgenie.db.tables import Asset, AssetGroup

    asset = session.exec(select(Asset).where(Asset.digest == link.asset_digest)).first()
    asset_group = (
        session.exec(select(AssetGroup).where(AssetGroup.id == asset.asset_group_id)).first()
        if asset
        else None
    )
    if not asset_group:
        return None
    if link.mode == "archive":
        from refgenie.utils.staging import staged_archive_path

        return staged_archive_path(
            genome_stage_folder, asset_group.genome_digest, asset_group.name, asset.digest
        )
    return genome_stage_folder / asset_group.genome_digest / asset_group.name / asset.name


def _push_folder_sync(results, session, genome_stage_folder, dry_run) -> tuple[int, int, int]:
    """Push via folder-level sync, marking links pushed individually.

    Groups unpushed links by remote and executes one sync command per remote.
    A link is only marked pushed if its staged file actually exists locally --
    a sync that exits 0 after silently skipping a missing file must not mark
    that link pushed forever -- and, for S3 remotes, if the uploaded object
    verifies afterward (see ``_verify_s3_object``). A remote with no
    push_command is a skip, not a failure. Returns
    ``(pushed_count, skipped_count, failed_count)``.
    """
    from collections import defaultdict

    by_remote = defaultdict(list)
    for link, remote, staged_asset in results:
        by_remote[remote.id].append((link, remote, staged_asset))

    pushed = 0
    skipped = 0
    failed = 0
    for items in by_remote.values():
        remote = items[0][1]

        if remote.push_command is None:
            logger.warning(f"Remote '{remote.description}' has no push_command. Skipping.")
            skipped += len(items)
            continue

        sync_cmd = remote.push_command
        sync_cmd = sync_cmd.replace("{genome_stage_folder}", str(genome_stage_folder))
        sync_cmd = sync_cmd.replace("{prefix}", remote.prefix)

        if dry_run:
            logger.info(f"[dry-run] Would execute: {sync_cmd}")
            logger.info(f"  Would mark {len(items)} assets as pushed to {remote.description}")
            pushed += len(items)
            continue

        logger.info(f"Syncing genome_stage_folder to {remote.description}: {sync_cmd}")
        try:
            result = subprocess.run(
                sync_cmd,
                shell=True,
                check=True,
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as e:
            logger.error(f"Folder sync failed for {remote.description}: {e}")
            if e.stderr:
                logger.error(f"stderr: {e.stderr.strip()}")
            failed += len(items)
            continue

        if result.stdout.strip():
            logger.debug(f"Sync stdout: {result.stdout.strip()}")

        s3_parts = _parse_s3_prefix(remote.prefix)

        marked = 0
        for link, _, staged_asset in items:
            local_path = _resolve_local_path(session, link, genome_stage_folder)
            if local_path is None or not local_path.exists():
                logger.error(
                    f"Not marking {link.asset_digest} mode={link.mode} as pushed: "
                    f"staged file missing locally ({local_path}), so the sync did not upload it."
                )
                failed += 1
                continue

            if s3_parts:
                bucket, base_path = s3_parts
                relative_path = local_path.relative_to(genome_stage_folder)
                s3_key = f"{base_path}/{relative_path}" if base_path else str(relative_path)
                if not _verify_s3_object(bucket, s3_key, staged_asset.tarball_size):
                    logger.error(
                        f"Not marking {link.asset_digest} mode={link.mode} as pushed: "
                        f"object did not verify at s3://{bucket}/{s3_key} after sync."
                    )
                    failed += 1
                    continue

            link.pushed = True
            session.add(link)
            marked += 1
        session.commit()
        pushed += marked
        logger.info(
            f"Synced and marked {marked}/{len(items)} assets as pushed to {remote.description}"
        )

    return pushed, skipped, failed


def _s3_head_object(bucket: str, key: str) -> dict | None:
    """Run ``aws s3api head-object`` and return the parsed response, or None.

    Returns the full JSON response (which carries ``ContentLength``, among
    other fields) rather than discarding it, so callers can verify an upload
    without a second network round-trip. Returns None if the object is
    absent, the call fails, or the aws CLI is not installed.
    """
    try:
        result = subprocess.run(
            ["aws", "s3api", "head-object", "--bucket", bucket, "--key", key],
            check=True,
            capture_output=True,
            text=True,
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None
    try:
        return json.loads(result.stdout)
    except json.JSONDecodeError:
        return None


def _verify_s3_object(bucket: str, key: str, expected_size: int | None) -> bool:
    """Confirm an object landed at s3://bucket/key after a push.

    Compares the remote's reported ``ContentLength`` against
    ``expected_size`` (``StagedAsset.tarball_size``). ``expected_size`` is
    None for mode="file" staged assets (only archive mode is a tarball with
    a recorded size); in that case this falls back to a presence-only check
    and logs at DEBUG that size could not be verified -- "no expected size"
    must never be silently treated as "verified".
    """
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


def _parse_s3_prefix(prefix: str) -> tuple[str, str] | None:
    """Parse an S3 prefix like 's3://bucket/path' into (bucket, path).

    Returns None if the prefix is not an S3 URL.
    """
    if not prefix.startswith("s3://"):
        return None
    rest = prefix[5:]
    if "/" in rest:
        bucket, path = rest.split("/", 1)
        return bucket, path.rstrip("/")
    return rest, ""


def _finish_push(
    pushed_count: int, skipped_count: int, failed_count: int, total: int, dry_run: bool
) -> None:
    """Log the push summary (omitting zero categories) and apply the exit contract.

    Shared by both push strategies so they agree on what counts as a skip vs
    a failure. A failed upload is always fatal. So is "nothing pushed, but
    something was skipped": the user asked to push and nothing landed, and
    exiting 0 on that is the worst outcome.
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
    """Push staged assets to cloud remotes.

    Queries RemoteAssetLink WHERE pushed=False, uploads each asset using
    the Remote's push_command template, then sets pushed=True.

    For S3-typed remotes with content-addressed archive mode, checks if the
    object already exists and skips the upload (idempotent). On dry-run,
    logs what would be skipped without any DB writes. For S3-typed remotes,
    a freshly uploaded object is verified (presence, and size when the
    staged record carries one) before its link is marked pushed.

    A remote with no push_command is a SKIP. A failed upload or a failed
    post-upload verification (nonzero exit, missing DB records, wrong
    object size) is a FAILURE. The command exits non-zero if any upload
    failed, or if nothing pushed while something was skipped.
    """
    from sqlmodel import select

    from refgenie.db.tables import (
        Asset,
        AssetGroup,
        Configuration,
        Remote,
        RemoteAssetLink,
        StagedAsset,
    )

    with refgenie._database_session as session:
        # Build query for unpushed links
        query = (
            select(RemoteAssetLink, Remote, StagedAsset)
            .join(Remote, RemoteAssetLink.remote_id == Remote.id)
            .join(
                StagedAsset,
                (StagedAsset.asset_digest == RemoteAssetLink.asset_digest)
                & (StagedAsset.mode == RemoteAssetLink.mode),
            )
            .where(RemoteAssetLink.pushed == False)  # noqa: E712
        )

        # Optional filters
        if cmd.remote:
            try:
                remote_id = int(cmd.remote)
                query = query.where(Remote.id == remote_id)
            except ValueError:
                query = query.where(Remote.description == cmd.remote)

        if cmd.genome:
            genome_digest = refgenie.alias.resolve(cmd.genome)
            query = (
                query.join(Asset, Asset.digest == RemoteAssetLink.asset_digest)
                .join(AssetGroup, AssetGroup.id == Asset.asset_group_id)
                .where(AssetGroup.genome_digest == genome_digest)
            )

        results = session.exec(query).all()

        if not results:
            logger.info("Nothing to push. All remote-asset links are up to date.")
            return

        # Get genome_stage_folder for path resolution
        config = session.exec(select(Configuration)).first()
        if not config or not config.genome_stage_folder:
            raise ValueError("genome_stage_folder is not configured")
        genome_stage_folder = Path(config.genome_stage_folder)

        # Handle folder_sync strategy
        if cmd.strategy == "folder_sync":
            pushed_count, skipped_count, failed_count = _push_folder_sync(
                results, session, genome_stage_folder, cmd.dry_run
            )
            _finish_push(pushed_count, skipped_count, failed_count, len(results), cmd.dry_run)
            return

        # Per-asset push strategy (default)
        pushed_count = 0
        skipped_count = 0
        failed_count = 0

        for link, remote, staged_asset in results:
            if remote.push_command is None:
                logger.warning(
                    f"Remote '{remote.description}' (id={remote.id}) has no push_command. "
                    f"Skipping {link.asset_digest} mode={link.mode}."
                )
                skipped_count += 1
                continue

            local_path = _resolve_local_path(session, link, genome_stage_folder)
            if local_path is None:
                logger.error(f"Asset records for {link.asset_digest} not found in DB. Cannot push.")
                failed_count += 1
                continue

            relative_path = local_path.relative_to(genome_stage_folder)

            # bucket/s3_key are also reused below to verify the upload landed.
            s3_parts = _parse_s3_prefix(remote.prefix)
            bucket, s3_key = (None, None)
            if s3_parts:
                bucket, base_path = s3_parts
                s3_key = f"{base_path}/{relative_path}" if base_path else str(relative_path)

            # Skip-if-present: for S3 remotes with archive mode, check if the content-addressed
            # object already exists. This enables idempotent re-pushes and free renames.
            if link.mode == "archive" and s3_key is not None:
                if _s3_head_object(bucket, s3_key) is not None:
                    if cmd.dry_run:
                        logger.info(
                            f"[dry-run] Would skip (already present): s3://{bucket}/{s3_key}"
                        )
                    else:
                        link.pushed = True
                        session.add(link)
                        session.commit()
                        logger.info(
                            f"Skipped (already present): s3://{bucket}/{s3_key} "
                            f"[{link.asset_digest}]"
                        )
                    pushed_count += 1
                    continue

            push_cmd = remote.push_command
            push_cmd = push_cmd.replace("{local_path}", str(local_path))
            push_cmd = push_cmd.replace("{relative_path}", str(relative_path))
            push_cmd = push_cmd.replace("{prefix}", remote.prefix)

            if cmd.dry_run:
                logger.info(f"[dry-run] Would execute: {push_cmd}")
                logger.info(
                    f"  asset_digest={link.asset_digest} mode={link.mode} -> {remote.description}"
                )
                pushed_count += 1
                continue

            logger.info(f"Pushing {local_path} to {remote.description}: {push_cmd}")
            try:
                result = subprocess.run(
                    push_cmd,
                    shell=True,
                    check=True,
                    capture_output=True,
                    text=True,
                )
                if result.stdout.strip():
                    logger.debug(f"Push stdout: {result.stdout.strip()}")

                if s3_key is not None and not _verify_s3_object(
                    bucket, s3_key, staged_asset.tarball_size
                ):
                    logger.error(
                        f"Not marking {link.asset_digest} mode={link.mode} as pushed: "
                        f"object did not verify at s3://{bucket}/{s3_key} after upload."
                    )
                    failed_count += 1
                    continue

                link.pushed = True
                session.add(link)
                session.commit()
                pushed_count += 1
                logger.info(f"Pushed {link.asset_digest} mode={link.mode} to {remote.description}")

            except subprocess.CalledProcessError as e:
                logger.error(f"Push failed for {link.asset_digest} to {remote.description}: {e}")
                if e.stderr:
                    logger.error(f"stderr: {e.stderr.strip()}")
                failed_count += 1

        _finish_push(pushed_count, skipped_count, failed_count, len(results), cmd.dry_run)
