"""
Staging, reached as ``rgc.stage``: making a built asset servable.

Staging writes an asset's serving form into the genome stage folder -- a
directory symlink for file mode, a tarball for archive mode -- and records a
``StagedAsset`` row. ``refgenie.utils.staging`` holds only the stage-folder path
convention; ``refgenie stage`` (``cli/commands/stage.py``) is the CLI over this.
"""

import os
import shutil
from pathlib import Path
from collections.abc import Iterable

from rich.table import Table
from sqlmodel import join, select
from sqlalchemy.engine import Engine
from sqlalchemy.orm import selectinload

from refgenie.db.tables import StagedAsset, Asset, AssetGroup
from refgenie.exceptions import MissingAssetError, MissingRemoteError, MissingStagedAssetError
from refgenie.logger import logger
from refgenie.managers.base import ResourceManager
from refgenie.managers.queries import one_or_raise
from refgenie.managers.remote import RemoteManager
from refgenie.utils.tables import build_table
from refgenie.utils.tarball import (
    get_external_symlinks,
    read_asset_dir_contents,
    read_build_commands,
    tar,
)
from refgenie.utils.digest import checksum
from refgenie.utils.staging import staged_archive_path
from ubiquerg import filesize_to_str as file_size


class StageManager(ResourceManager):
    """
    A manager for staged assets. Staging prepares assets for serving and/or
    cloud push by creating entries in genome_stage_folder.

    File-mode: creates a directory symlink from genome_stage_folder to genome_folder.
    Archive-mode: creates a .tgz tarball in genome_stage_folder.
    """

    def __init__(self, database_engine: Engine, remote_manager: RemoteManager):
        super().__init__(database_engine)
        self._remote = remote_manager

    def create(
        self,
        asset: Asset,
        genome_folder: Path,
        genome_stage_folder: Path | None,
        push_to: list[str] | None = None,
        build_dir: Path | None = None,
    ) -> list[StagedAsset]:
        """
        Stage an asset according to its serving modes.

        - If "archive" in modes: create .tgz tarball, create StagedAsset(mode="archive")
        - If "file" in modes: link the files, create StagedAsset(mode="file")
        - If both: do both (two StagedAsset records)
        - If "none" only: skip entirely, return empty list

        File-mode staging creates a directory of per-file symlinks under
        genome_stage_folder/.../asset_name into genome_folder. This avoids disk duplication. The OS follows symlinks transparently for serving,
        and cloud sync tools (aws s3 sync, etc.) follow symlinks by default.

        Args:
            asset: The asset to stage.
            genome_folder: Path to the genome folder.
            genome_stage_folder: Path to the staging folder.
            push_to: Optional remote names or ids to create push intent records for.
            build_dir: Optional path to the asset's build bookkeeping directory,
                used to record the build commands. Resolving it is the caller's
                job; pulled assets have none.

        Returns:
            The asset's StagedAsset records (0, 1, or 2 items).
        """
        logger.info(f"Staging asset: {asset}")

        if genome_stage_folder is None:
            raise ValueError("genome_stage_folder must be set to stage assets")

        if asset.digest is None:
            raise ValueError(f"Can't stage. Can't identify asset digest for: {asset}")

        modes = self._serving_modes(asset)
        if set(modes) == {"none"}:
            logger.info(f"Serving mode is 'none' for {asset}. Skipping staging.")
            return []

        # Existing records are kept; each mode's helper repairs what is stale
        existing = self._get_existing_staged(asset.digest)
        input_dir = genome_folder / asset.path
        build_commands = read_build_commands(build_dir) if build_dir else []

        results = list(existing)
        if "archive" in modes:
            staged = self._stage_archive(
                asset, input_dir, genome_stage_folder, existing, build_commands
            )
            if staged is not None:
                results = [r for r in results if r.id != staged.id]
                results.append(staged)
        if "file" in modes:
            staged = self._stage_files(
                asset, input_dir, genome_stage_folder, existing, build_commands
            )
            if staged is not None:
                results.append(staged)

        if push_to and results:
            self._record_push_intent(push_to, results)

        return results

    def _serving_modes(self, asset: Asset) -> list[str]:
        """
        Read an asset's serving modes, reloading it if it is detached.

        Args:
            asset: The asset, possibly detached from its session.

        Returns:
            list[str]: The asset's serving modes.
        """
        try:
            return asset.serving_modes
        except Exception:
            with self._database_session as session:
                loaded_asset = one_or_raise(
                    session,
                    select(Asset)
                    .options(selectinload(Asset.asset_group).selectinload(AssetGroup.asset_class))
                    .where(Asset.digest == asset.digest),
                    MissingAssetError(digest=asset.digest),
                    unique=True,
                )
                return loaded_asset.serving_modes

    def _stage_archive(
        self,
        asset: Asset,
        input_dir: Path,
        genome_stage_folder: Path,
        existing: list[StagedAsset],
        build_commands: list[str],
    ) -> StagedAsset | None:
        """
        Write the asset's tarball and its archive-mode record.

        The tarball is content-addressed by asset digest:
        ``{genome_digest}/{group}/{asset_digest}.tgz``. An existing record is a
        claim that the tarball is there and is that digest; it is not evidence.
        Staging was interrupted between writing the tarball and committing the
        row often enough to make the claim worth checking, and a row pointing at
        a missing or wrong tarball publishes a download that fails for every
        client. So an existing record is repaired in place when its tarball is
        not intact.

        The tarball is written to a ``.partial`` sibling first and only moved
        onto the final content-addressed path (via ``os.replace``, atomic on
        the same filesystem) once it has been fully written and checksummed.
        The final path therefore only ever holds a complete tarball: DRS and
        v4 byte-serving routes require a ``StagedAsset`` row before serving
        from this path (see ``require_published``), and that row is committed
        only after the swap, so a tarball with no row is never served, and a
        row never points at a partial file. In the repair case the previous
        (possibly stale) tarball stays at the final path, still servable,
        until the new one replaces it.

        Args:
            asset: The asset to stage.
            input_dir: The asset's content directory.
            genome_stage_folder: Path to the staging folder.
            existing: The asset's existing staged records.
            build_commands: The build commands to record.

        Returns:
            StagedAsset | None: The created or repaired record, or None when an
                intact archive was already staged.
        """
        target_file = staged_archive_path(
            genome_stage_folder,
            asset.asset_group.genome_digest,
            asset.asset_group.name,
            asset.digest,
        )
        existing_archive = next((sa for sa in existing if sa.mode == "archive"), None)
        if existing_archive is not None and self._archive_is_intact(target_file, existing_archive):
            return None

        target_file.parent.mkdir(parents=True, exist_ok=True)
        # Exclude colocation symlinks (external symlinks to parent assets) from
        # the archive; they are recreated after pull.
        exclude_files = get_external_symlinks(input_dir)
        directory_contents = read_asset_dir_contents(input_dir, exclude_files=set(exclude_files))

        # Write and checksum in a .partial sibling, then atomically swap it
        # onto the final path. Nothing at the final path is ever a partial or
        # unchecksummed tarball -- a crash between tar() and the swap leaves
        # only a .partial file, which no route serves and the next stage
        # overwrites.
        partial_file = target_file.with_name(target_file.name + ".partial")
        try:
            tar(input_dir, partial_file, exclude_files=exclude_files)
            # Digest the tarball here, not in a before_insert listener. A
            # multi-GB checksum inside the insert flush holds the write lock
            # for its whole duration, for a value that has nothing to do with
            # the transaction.
            tarball_digest = checksum(str(partial_file))
            tarball_size = partial_file.stat().st_size
        except BaseException:
            partial_file.unlink(missing_ok=True)
            raise
        os.replace(partial_file, target_file)

        with self._database_session as session:
            if existing_archive is not None:
                sa = session.get(StagedAsset, existing_archive.id)
                sa.build_commands = build_commands
                sa.directory_contents = directory_contents
                sa.tarball_digest = tarball_digest
                sa.tarball_size = tarball_size
                action = "Repaired"
            else:
                sa = StagedAsset(
                    asset_digest=asset.digest,
                    mode="archive",
                    build_commands=build_commands,
                    directory_contents=directory_contents,
                    tarball_digest=tarball_digest,
                    tarball_size=tarball_size,
                )
                action = "Created"
            session.add(sa)
            session.commit()
            session.refresh(sa)
        logger.info(f"{action} archive StagedAsset for {asset}")
        return sa

    def _stage_files(
        self,
        asset: Asset,
        input_dir: Path,
        genome_stage_folder: Path,
        existing: list[StagedAsset],
        build_commands: list[str],
    ) -> StagedAsset | None:
        """
        Link the asset's files into the stage folder and record file mode.

        The stage entry is a directory of per-file symlinks into the asset's
        real output files, excluding colocation symlinks (external symlinks to
        parent assets). This mirrors archive mode, which strips those from the
        tarball. A whole-directory symlink would expose them, and a relative
        colocation symlink (e.g. ../../fasta/default/x.fa) then fails to resolve
        once the asset dir is reached through the stage entry, which breaks
        file-mode push (`aws s3 sync --follow-symlinks` cannot find its target).
        Colocation symlinks stay in the source build dir for consumers; they are
        not part of the published or served asset.

        The links are reconciled on every call, even when the record exists,
        because they are derived from the asset and a killed run can leave a
        half-linked directory.

        Args:
            asset: The asset to stage.
            input_dir: The asset's content directory.
            genome_stage_folder: Path to the staging folder.
            existing: The asset's existing staged records.
            build_commands: The build commands to record.

        Returns:
            StagedAsset | None: The created record, or None when file mode was
                already staged.
        """
        link_dir = (
            genome_stage_folder
            / asset.asset_group.genome_digest
            / asset.asset_group.name
            / asset.name
        )
        link_dir.parent.mkdir(parents=True, exist_ok=True)
        external = set(get_external_symlinks(input_dir))
        self._reconcile_stage_links(link_dir, input_dir, external)
        directory_contents = read_asset_dir_contents(input_dir, exclude_files=external)

        if any(sa.mode == "file" for sa in existing):
            return None
        with self._database_session as session:
            sa = StagedAsset(
                asset_digest=asset.digest,
                mode="file",
                build_commands=build_commands,
                directory_contents=directory_contents,
                # tarball_digest and tarball_size stay None for file mode
            )
            session.add(sa)
            session.commit()
            session.refresh(sa)
        logger.info(f"Created file-serving StagedAsset for {asset}")
        return sa

    def _record_push_intent(self, push_to: list[str], staged: list[StagedAsset]) -> None:
        """
        Link each staged record to each requested remote.

        Idempotent: a persistent catalog re-stages already-staged assets every
        nightly run, so a link for this (remote, asset, mode) may already exist.
        ``exist_ok`` keeps it as is, preserving its pushed state (an
        already-pushed asset must stay pushed, not reset). An unknown remote is
        skipped with a warning.

        Args:
            push_to: Remote names or ids.
            staged: The staged records to push.
        """
        for ref in push_to:
            created = 0
            try:
                for sa in staged:
                    _, new = self._remote.link(ref, sa.asset_digest, sa.mode, exist_ok=True)
                    created += new
            except MissingRemoteError:
                logger.warning(f"Remote '{ref}' not found. Skipping push intent.")
                continue
            logger.info(
                f"Push intent for remote '{ref}': "
                f"{created} new, {len(staged) - created} already tracked"
            )

    @staticmethod
    def _archive_is_intact(target_file: Path, staged: StagedAsset) -> bool:
        """
        Whether the tarball on disk is the one ``staged`` describes.

        Size is checked before the digest because a truncated tarball -- what a
        killed ``tar`` leaves behind -- is decidable from a ``stat``, and because
        a mismatch there makes the read pointless.

        Args:
            target_file: The content-addressed tarball path.
            staged: The existing archive-mode record.

        Returns:
            bool: True when the tarball exists and matches the record.
        """
        if not target_file.is_file():
            logger.warning(f"Staged tarball is missing; re-creating it: {target_file}")
            return False
        if target_file.stat().st_size != staged.tarball_size:
            logger.warning(f"Staged tarball size differs from the record: {target_file}")
            return False
        if checksum(str(target_file)) != staged.tarball_digest:
            logger.warning(f"Staged tarball digest differs from the record: {target_file}")
            return False
        return True

    @staticmethod
    def _reconcile_stage_links(link_dir: Path, input_dir: Path, external: set[str]) -> None:
        """
        Make ``link_dir`` hold exactly one symlink per served file in ``input_dir``.

        The served set is derived from the asset, so this is idempotent and safe
        to run on a directory in any state: absent, complete, or half-written by
        an interrupted stage.

        Args:
            link_dir: The staging directory of per-file symlinks.
            input_dir: The asset's content directory.
            external: Names of colocation symlinks, which are not served.
        """
        if link_dir.is_symlink():
            # Never walk through a symlinked stage dir: pruning stale entries
            # below would then delete the asset's own files.
            link_dir.unlink()
        link_dir.mkdir(parents=True, exist_ok=True)

        expected = {
            child.name: child
            for child in sorted(input_dir.iterdir(), key=lambda p: p.name)
            if child.name not in external
        }
        repaired = 0
        for entry in link_dir.iterdir():
            wanted = expected.get(entry.name)
            if wanted is not None and entry.is_symlink() and Path(os.readlink(entry)) == wanted:
                continue
            # Stale, wrongly-targeted, or not a symlink at all.
            if entry.is_symlink() or entry.is_file():
                entry.unlink()
            else:
                shutil.rmtree(entry)
            repaired += 1
        for name, child in expected.items():
            if not (link_dir / name).is_symlink():
                (link_dir / name).symlink_to(child)
                repaired += 1

        if repaired:
            logger.info(
                f"Reconciled file-serving stage dir: {link_dir} "
                f"({len(expected)} link(s), {repaired} written or removed; "
                f"excluded {len(external)} colocation symlink(s))"
            )

    def remove(self, asset_digest: str):
        """
        Remove all StagedAsset records for an asset (unstage).
        The before_delete event handler handles file/symlink cleanup.

        Args:
            asset_digest: The digest of the asset.
        """
        with self._database_session as session:
            records = session.exec(
                select(StagedAsset).where(StagedAsset.asset_digest == asset_digest)
            ).all()
            if not records:
                logger.warning(f"No StagedAsset records found for {asset_digest}")
                return
            for record in records:
                session.delete(record)
            session.commit()
        logger.info(f"Removed StagedAsset records for {asset_digest}")

    def get_by_asset_digest(self, asset_digest: str) -> list[StagedAsset]:
        """
        Get all staged asset records by asset digest.

        Args:
            asset_digest: The digest of the asset.

        Returns:
            list[StagedAsset]: The staged asset records.
        """
        with self._database_session as session:
            results = list(
                session.exec(
                    select(StagedAsset)
                    .options(
                        selectinload(StagedAsset.asset)
                        .selectinload(Asset.asset_group)
                        .selectinload(AssetGroup.genome)
                    )
                    .where(StagedAsset.asset_digest == asset_digest)
                ).all()
            )
            if not results:
                raise MissingStagedAssetError(f"{asset_digest=}")
            return results

    def exists_by_asset_digest(self, asset_digest: str) -> bool:
        """
        Check if any staged asset records exist for the given asset digest.

        Args:
            asset_digest: The digest of the asset.

        Returns:
            bool: True if staged asset records exist, False otherwise.
        """
        with self._database_session as session:
            return bool(
                session.exec(
                    select(StagedAsset).where(
                        StagedAsset.asset_digest == asset_digest,
                    )
                ).first()
            )

    def increment_download_count(self, asset_digest: str, mode: str, count: int = 1):
        """
        Increment the download count of a staged asset.

        Args:
            asset_digest: The digest of the asset.
            mode: The mode ("file" or "archive").
            count: The count to increment by.
        """
        with self._database_session as session:
            sa = one_or_raise(
                session,
                select(StagedAsset).where(
                    StagedAsset.asset_digest == asset_digest,
                    StagedAsset.mode == mode,
                ),
                MissingStagedAssetError(f"{asset_digest=} {mode=}"),
            )
            sa.download_count += count
            session.commit()
        logger.info(
            f"Updated staged asset download count: {asset_digest=} {mode=} {sa.download_count=}"
        )

    def list_all(self) -> Iterable[tuple[StagedAsset, Asset]]:
        """
        List all staged assets and their associated assets.
        """

        with self._database_session as session:
            return (
                session.exec(select(StagedAsset, Asset).select_from(join(StagedAsset, Asset)))
                .unique()
                .all()
            )

    def table(self) -> Table:
        """
        Create a table of all staged assets.
        """
        rows = []
        with self._database_session as session:
            for sa in session.exec(select(StagedAsset)).all():
                asset = session.exec(select(Asset).where(Asset.digest == sa.asset_digest)).first()
                if asset is None:
                    continue
                rows.append(
                    (
                        sa.asset_digest,
                        f"{asset.asset_group.genome.digest}/{asset.asset_group.name}:{asset.name}",
                        sa.mode,
                        file_size(sa.tarball_size) if sa.tarball_size else "-",
                    )
                )

        return build_table("Staged Assets", ["Asset Digest", "Asset name", "Mode", "Size"], rows)

    def _get_existing_staged(self, asset_digest: str) -> list[StagedAsset]:
        with self._database_session as session:
            return list(
                session.exec(
                    select(StagedAsset).where(StagedAsset.asset_digest == asset_digest)
                ).all()
            )
