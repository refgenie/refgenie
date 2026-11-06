"""Queries on push remotes and the staged assets linked to each."""

from sqlalchemy.orm import selectinload
from sqlmodel import Session, select

from refgenie.db.queries.common import one_or_raise
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

#: Links have no id; their key is (remote, asset digest, mode).
_LINK_ORDER = (RemoteAssetLink.remote_id, RemoteAssetLink.asset_digest, RemoteAssetLink.mode)


def get(session: Session, ref: str | int) -> Remote:
    """
    The remote ``ref`` names.

    The one naming rule: an ``int``, or a string of digits, is the remote's id;
    any other string is its name.

    Raises:
        MissingRemoteError: If no remote matches.
    """
    if isinstance(ref, int) or ref.isdigit():
        statement = select(Remote).where(Remote.id == int(ref))
    else:
        statement = select(Remote).where(Remote.name == ref)
    return one_or_raise(session, statement, MissingRemoteError(ref), unique=True)


def find_by_name(session: Session, name: str) -> Remote | None:
    """The remote with this name, or None."""
    return session.exec(select(Remote).where(Remote.name == name)).unique().one_or_none()


def list_all(session: Session) -> list[Remote]:
    """Every remote, by id."""
    return list(session.exec(select(Remote).order_by(Remote.id)).unique().all())


def find_link(
    session: Session, remote_id: int, asset_digest: str, mode: str
) -> RemoteAssetLink | None:
    """The link of ``asset_digest`` in ``mode`` to a remote, or None."""
    return session.exec(
        select(RemoteAssetLink).where(
            RemoteAssetLink.remote_id == remote_id,
            RemoteAssetLink.asset_digest == asset_digest,
            RemoteAssetLink.mode == mode,
        )
    ).one_or_none()


def list_links(session: Session, remote_ids: list[int] | None = None) -> list[RemoteAssetLink]:
    """Every link, or the links of these remotes, by remote, asset digest and mode."""
    statement = select(RemoteAssetLink)
    if remote_ids is not None:
        statement = statement.where(RemoteAssetLink.remote_id.in_(remote_ids))
    return list(session.exec(statement.order_by(*_LINK_ORDER)).all())


def list_unpushed(
    session: Session, remote_id: int | None = None, genome_digest: str | None = None
) -> list[tuple[RemoteAssetLink, Remote, StagedAsset]]:
    """
    The links still waiting to be pushed, with their remote and staged asset.

    Each ``StagedAsset`` comes with its ``asset`` and ``asset.asset_group``
    loaded, so callers can derive staged paths after the session closes.
    """
    statement = (
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
    if remote_id is not None:
        statement = statement.where(Remote.id == remote_id)
    if genome_digest is not None:
        statement = (
            statement.join(Asset, Asset.digest == RemoteAssetLink.asset_digest)
            .join(AssetGroup, AssetGroup.id == Asset.asset_group_id)
            .where(AssetGroup.genome_digest == genome_digest)
        )
    statement = statement.order_by(*_LINK_ORDER)
    return [tuple(row) for row in session.exec(statement).all()]


def list_pushed_http(
    session: Session, asset_digest: str, modes: list[str]
) -> list[tuple[str, Remote, str]]:
    """
    ``(mode, remote, genome_stage_folder)`` for every http(s) remote this asset was pushed to.

    https first, then newest configuration first, so the first row for a mode
    is the preferred download.
    """
    return list(
        session.exec(
            select(RemoteAssetLink.mode, Remote, Configuration.genome_stage_folder)
            .join(Remote, RemoteAssetLink.remote_id == Remote.id)
            .join(Configuration, Remote.configuration_id == Configuration.id)
            .where(RemoteAssetLink.asset_digest == asset_digest)
            .where(RemoteAssetLink.mode.in_(modes))
            .where(RemoteAssetLink.pushed == True)  # noqa: E712
            .where(Remote.type.in_([RemoteType.https, RemoteType.http]))
            .order_by(Remote.type.desc(), Configuration.id.desc(), Remote.id)
        ).all()
    )
