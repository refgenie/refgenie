"""
Push targets and the staged assets queued for each.

"Remote" in refgenie means a place ``refgenie push`` uploads to (the ``Remote``
table). The servers a client pulls from are ``rgc.servers``, not this.
"""

import pathlib
from collections import defaultdict
from collections.abc import Mapping

from rich.table import Table
from sqlalchemy.orm import selectinload
from sqlmodel import Session, select

from refgenie.db.tables import (
    Asset,
    AssetGroup,
    Configuration,
    Remote,
    RemoteAssetLink,
    RemoteType,
    StagedAsset,
)
from refgenie.exceptions import MissingRemoteError
from refgenie.logger import logger
from refgenie.managers.base import ResourceManager
from refgenie.managers.configuration import latest_configuration
from refgenie.managers.queries import one_or_raise
from refgenie.models import GenomeDigest
from refgenie.utils.tables import build_table


def _remote_in_session(session: Session, ref: str | int) -> Remote:
    """
    The remote ``ref`` names, read through a caller's session.

    The one naming rule: an ``int``, or a string of digits, is the remote's id;
    any other string is its name.
    """
    if isinstance(ref, int) or ref.isdigit():
        query = select(Remote).where(Remote.id == int(ref))
    else:
        query = select(Remote).where(Remote.name == ref)
    return one_or_raise(session, query, MissingRemoteError(ref), unique=True)


def _check_name(name: str) -> None:
    """A name must be non-empty and not all digits, since digits read as an id."""
    if not name or name.isdigit():
        raise ValueError(f"Invalid remote name '{name}': it must not be empty or all digits.")


def _link_query(remote_id: int, asset_digest: str, mode: str):
    return select(RemoteAssetLink).where(
        RemoteAssetLink.remote_id == remote_id,
        RemoteAssetLink.asset_digest == asset_digest,
        RemoteAssetLink.mode == mode,
    )


class RemoteManager(ResourceManager):
    """
    Push targets (``Remote``) and which staged assets go to each (``RemoteAssetLink``).

    Every method that takes a remote ``ref`` accepts its id (an ``int`` or a
    string of digits) or its name.
    """

    # --- remotes ---

    def get(self, ref: str | int) -> Remote:
        """
        Get a remote by id or name.

        Args:
            ref: The remote's id (int or digit string) or name.

        Returns:
            Remote: The remote.

        Raises:
            MissingRemoteError: If no remote matches.
        """
        with self._database_session as session:
            return _remote_in_session(session, ref)

    def exists(self, ref: str | int) -> bool:
        """Whether a remote with this id or name exists."""
        try:
            self.get(ref)
        except MissingRemoteError:
            return False
        return True

    def list_all(self) -> list[Remote]:
        """All remotes, in id order."""
        with self._database_session as session:
            return list(session.exec(select(Remote).order_by(Remote.id)).unique().all())

    def add(
        self,
        name: str,
        type: RemoteType,
        prefix: str,
        push_command: str | None = None,
    ) -> Remote:
        """
        Add a remote. Multiple remotes of one type are allowed; names are unique.

        Args:
            name: The remote's name. Must not be all digits, since that reads as an id.
            type: The remote type.
            prefix: The remote prefix (bucket or base URL).
            push_command: Optional shell command template for pushing assets.

        Returns:
            Remote: The added remote.

        Raises:
            ValueError: If the name is empty, all digits, or already taken.
        """
        _check_name(name)
        with self._database_session as session:
            if session.exec(select(Remote).where(Remote.name == name)).first():
                raise ValueError(f"A remote named '{name}' already exists.")
            remote = Remote(
                type=RemoteType(type),
                prefix=prefix,
                name=name,
                push_command=push_command,
                configuration_id=latest_configuration(session).id,
            )
            session.add(remote)
            session.commit()
            session.refresh(remote)
            logger.info(f"Added remote '{name}' (id={remote.id})")
            return remote

    def upsert(
        self,
        name: str,
        type: RemoteType,
        prefix: str,
        push_command: str | None = None,
    ) -> Remote:
        """
        Update the remote with this name, or add it.

        Args:
            name: The remote's name.
            type: The remote type.
            prefix: The remote prefix.
            push_command: Optional shell command template for pushing assets.

        Returns:
            Remote: The updated or added remote.
        """
        _check_name(name)
        with self._database_session as session:
            remote = session.exec(select(Remote).where(Remote.name == name)).first()
            created = remote is None
            if created:
                remote = Remote(name=name, configuration_id=latest_configuration(session).id)
            remote.type = RemoteType(type)
            remote.prefix = prefix
            remote.push_command = push_command
            session.add(remote)
            session.commit()
            session.refresh(remote)
            logger.info(f"{'Added' if created else 'Updated'} remote '{name}' (id={remote.id})")
            return remote

    def remove(self, ref: str | int) -> None:
        """
        Remove a remote by id or name.

        Raises:
            MissingRemoteError: If no remote matches.
        """
        with self._database_session as session:
            remote = _remote_in_session(session, ref)
            session.delete(remote)
            session.commit()
            logger.info(f"Removed remote '{remote.name}' (id={remote.id})")

    def table(self) -> Table:
        """The remotes, with pushed and unpushed link counts."""
        with self._database_session as session:
            remotes = session.exec(select(Remote).order_by(Remote.id)).unique().all()
            links = session.exec(select(RemoteAssetLink)).all()
        pushed_counts = defaultdict(int)
        unpushed_counts = defaultdict(int)
        for link in links:
            if link.pushed:
                pushed_counts[link.remote_id] += 1
            else:
                unpushed_counts[link.remote_id] += 1

        return build_table(
            "Remotes",
            ["ID", "Name", "Type", "Prefix", "Push Command", "Pushed", "Unpushed"],
            [
                (
                    str(remote.id),
                    remote.name,
                    remote.type.value,
                    remote.prefix,
                    remote.push_command or "-",
                    str(pushed_counts.get(remote.id, 0)),
                    str(unpushed_counts.get(remote.id, 0)),
                )
                for remote in remotes
            ],
        )

    def status(self, ref: str | int | None = None) -> dict[int, dict]:
        """
        Push status per remote.

        Returns a dict of remote_id ->
        {"remote": Remote, "pushed": [links], "unpushed": [links]}.
        Remotes with no asset links are included, with empty lists.

        Args:
            ref: Optional remote id or name to report on alone.

        Raises:
            MissingRemoteError: If ``ref`` names no remote.
        """
        with self._database_session as session:
            if ref is None:
                remotes = session.exec(select(Remote).order_by(Remote.id)).unique().all()
            else:
                remotes = [_remote_in_session(session, ref)]
            if not remotes:
                return {}
            by_remote = {r.id: {"remote": r, "pushed": [], "unpushed": []} for r in remotes}
            links = session.exec(
                select(RemoteAssetLink).where(RemoteAssetLink.remote_id.in_(by_remote))
            ).all()
        for link in links:
            by_remote[link.remote_id]["pushed" if link.pushed else "unpushed"].append(link)
        return by_remote

    # --- links ---

    def link(
        self,
        ref: str | int,
        asset_digest: str,
        mode: str,
        *,
        pushed: bool = False,
        exist_ok: bool = False,
    ) -> tuple[RemoteAssetLink, bool]:
        """
        Record that a staged asset goes to a remote.

        Args:
            ref: The remote's id or name.
            asset_digest: The asset digest.
            mode: The staging mode ("archive" or "file").
            pushed: Whether the asset is already pushed (default False: an intent).
            exist_ok: Return an existing link unchanged instead of raising. Its
                ``pushed`` state is never reset: re-staging an already pushed
                asset must leave it pushed.

        Returns:
            (link, created): The link, and whether this call created it.

        Raises:
            MissingRemoteError: If ``ref`` names no remote.
            ValueError: If the (asset_digest, mode) StagedAsset does not exist,
                or the link exists and ``exist_ok`` is False.
        """
        with self._database_session as session:
            remote = _remote_in_session(session, ref)
            one_or_raise(
                session,
                select(StagedAsset).where(
                    StagedAsset.asset_digest == asset_digest,
                    StagedAsset.mode == mode,
                ),
                ValueError(
                    f"StagedAsset not found for {asset_digest=}, {mode=}. "
                    f"Asset must be staged before linking to a remote."
                ),
            )
            existing = session.exec(_link_query(remote.id, asset_digest, mode)).one_or_none()
            if existing is not None:
                if exist_ok:
                    return existing, False
                raise ValueError(
                    f"{asset_digest} ({mode}) is already linked to remote '{remote.name}'"
                )
            link = RemoteAssetLink(
                remote_id=remote.id, asset_digest=asset_digest, mode=mode, pushed=pushed
            )
            session.add(link)
            session.commit()
            session.refresh(link)
            logger.info(
                f"Linked {asset_digest} ({mode}) to remote '{remote.name}' (pushed={pushed})"
            )
            return link, True

    def unlink(self, ref: str | int, asset_digest: str, mode: str) -> None:
        """
        Remove the link between a staged asset and a remote.

        Raises:
            MissingRemoteError: If ``ref`` names no remote.
            ValueError: If there is no such link.
        """
        with self._database_session as session:
            remote = _remote_in_session(session, ref)
            link = one_or_raise(
                session,
                _link_query(remote.id, asset_digest, mode),
                ValueError(f"No link found for remote '{ref}', {asset_digest=}, {mode=}"),
            )
            session.delete(link)
            session.commit()
            logger.info(f"Unlinked {asset_digest} ({mode}) from remote '{remote.name}'")

    def mark_pushed(self, ref: str | int, asset_digest: str, mode: str) -> None:
        """
        Mark a link pushed. ``refgenie push`` calls this after an upload succeeds.

        Raises:
            MissingRemoteError: If ``ref`` names no remote.
            ValueError: If there is no such link.
        """
        with self._database_session as session:
            remote = _remote_in_session(session, ref)
            link = one_or_raise(
                session,
                _link_query(remote.id, asset_digest, mode),
                ValueError(f"No link found for remote '{ref}', {asset_digest=}, {mode=}"),
            )
            link.pushed = True
            session.add(link)
            session.commit()
            logger.info(f"Marked pushed: {asset_digest} ({mode}) on remote '{remote.name}'")

    def unpushed(
        self,
        ref: str | int | None = None,
        *,
        genome_digest: GenomeDigest | None = None,
    ) -> list[tuple[RemoteAssetLink, Remote, StagedAsset]]:
        """
        The links still waiting to be pushed, with their remote and staged asset.

        Each ``StagedAsset`` comes with its ``asset`` and ``asset.asset_group``
        loaded, so callers can derive staged paths after the session closes.

        Args:
            ref: Optional remote id or name to limit to.
            genome_digest: Optional genome to limit to.

        Raises:
            MissingRemoteError: If ``ref`` names no remote.
        """
        with self._database_session as session:
            query = (
                select(RemoteAssetLink, Remote, StagedAsset)
                .join(Remote, RemoteAssetLink.remote_id == Remote.id)
                .join(
                    StagedAsset,
                    (StagedAsset.asset_digest == RemoteAssetLink.asset_digest)
                    & (StagedAsset.mode == RemoteAssetLink.mode),
                )
                .where(RemoteAssetLink.pushed == False)  # noqa: E712
                .options(selectinload(StagedAsset.asset).selectinload(Asset.asset_group))
            )
            if ref is not None:
                query = query.where(Remote.id == _remote_in_session(session, ref).id)
            if genome_digest is not None:
                query = (
                    query.join(Asset, Asset.digest == RemoteAssetLink.asset_digest)
                    .join(AssetGroup, AssetGroup.id == Asset.asset_group_id)
                    .where(AssetGroup.genome_digest == genome_digest)
                )
            return [tuple(row) for row in session.exec(query).all()]

    def pushed_urls(
        self, asset_digest: str, local_paths: Mapping[str, pathlib.Path]
    ) -> list[tuple[Remote, str, str]]:
        """
        Where this asset can be downloaded from, on remotes it was pushed to.

        Only http and https remotes give a URL. The URL is the remote prefix
        plus the local path relative to the stage folder. Rows come https
        first, then newest configuration first, so the first row for a mode is
        the preferred download.

        Args:
            asset_digest: The asset digest.
            local_paths: mode -> the local staged path to translate (the
                tarball for "archive"; the asset directory, or a file in it, for
                "file"). Only these modes are considered.

        Returns:
            (remote, mode, url) for every pushed link that yields a URL.
        """
        if not local_paths:
            return []
        with self._database_session as session:
            rows = session.exec(
                select(RemoteAssetLink.mode, Remote, Configuration.genome_stage_folder)
                .join(Remote, RemoteAssetLink.remote_id == Remote.id)
                .join(Configuration, Remote.configuration_id == Configuration.id)
                .where(RemoteAssetLink.asset_digest == asset_digest)
                .where(RemoteAssetLink.mode.in_(list(local_paths)))
                .where(RemoteAssetLink.pushed == True)  # noqa: E712
                .where(Remote.type.in_([RemoteType.https, RemoteType.http]))
                .order_by(Remote.type.desc(), Configuration.id.desc())
            ).all()

        urls = []
        for mode, remote, genome_stage_folder in rows:
            local_path = local_paths[mode]
            try:
                relative_path = str(local_path.relative_to(pathlib.Path(genome_stage_folder)))
            except ValueError:
                logger.warning(
                    f"Cannot determine remote path for {asset_digest=}. "
                    f"{local_path=} is not relative to {genome_stage_folder=}."
                )
                continue
            urls.append((remote, mode, f"{remote.prefix.rstrip('/')}/{relative_path.lstrip('/')}"))
        return urls
