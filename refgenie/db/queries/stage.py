"""Queries on staged assets: the publication records of the serving modes."""

from sqlalchemy.orm import selectinload
from sqlmodel import Session, select

from refgenie.db.queries.common import one_or_raise
from refgenie.db.queries.paging import paginate
from refgenie.db.queries.search import search_condition
from refgenie.db.tables import Asset, AssetGroup, StagedAsset
from refgenie.exceptions import MissingStagedAssetError
from refgenie.utils.search import Search

#: Fields a staged asset search may match.
SEARCH_FIELDS = ["asset_digest", "mode"]

_ORDER = (StagedAsset.asset_digest, StagedAsset.mode, StagedAsset.id)

#: The staged row's asset with its group, genome and asset class: enough for
#: ``registry_path`` and ``serving_modes``.
_WITH_ASSET = (
    selectinload(StagedAsset.asset).selectinload(Asset.asset_group).selectinload(AssetGroup.genome),
    selectinload(StagedAsset.asset)
    .selectinload(Asset.asset_group)
    .selectinload(AssetGroup.asset_class),
)


def get_by_id(session: Session, id: int) -> StagedAsset | None:
    """The staged asset row with this id, or None."""
    return session.get(StagedAsset, id)


def find(session: Session, asset_digest: str, mode: str) -> StagedAsset | None:
    """The staged row of ``asset_digest`` in ``mode``, or None if it is not published so."""
    statement = select(StagedAsset).where(
        StagedAsset.asset_digest == asset_digest, StagedAsset.mode == mode
    )
    return session.exec(statement).one_or_none()


def get(session: Session, asset_digest: str, mode: str) -> StagedAsset:
    """
    The staged row of ``asset_digest`` in ``mode``.

    Raises:
        MissingStagedAssetError: If there is none.
    """
    return one_or_raise(
        session,
        select(StagedAsset).where(
            StagedAsset.asset_digest == asset_digest, StagedAsset.mode == mode
        ),
        MissingStagedAssetError(f"{asset_digest=} {mode=}"),
    )


def list_for_asset(
    session: Session, asset_digest: str, *, with_asset: bool = False
) -> list[StagedAsset]:
    """Every staged row of one asset, by mode."""
    statement = select(StagedAsset).where(StagedAsset.asset_digest == asset_digest)
    if with_asset:
        statement = statement.options(*_WITH_ASSET)
    return list(session.exec(statement.order_by(*_ORDER)).all())


def exists_for_asset(session: Session, asset_digest: str) -> bool:
    """Whether the asset is staged in any mode."""
    statement = select(StagedAsset.id).where(StagedAsset.asset_digest == asset_digest)
    return session.exec(statement).first() is not None


def list_all(session: Session) -> list[StagedAsset]:
    """Every staged row with its asset loaded, by asset digest and mode."""
    return list(session.exec(select(StagedAsset).options(*_WITH_ASSET).order_by(*_ORDER)).all())


def list_archives(
    session: Session, *, asset_digest: str | None = None, genome_digest: str | None = None
) -> list[StagedAsset]:
    """Archive-mode staged rows with their asset loaded, optionally of one asset or genome."""
    statement = select(StagedAsset).options(*_WITH_ASSET).where(StagedAsset.mode == "archive")
    if asset_digest:
        statement = statement.where(StagedAsset.asset_digest == asset_digest)
    if genome_digest:
        statement = (
            statement.join(Asset, Asset.digest == StagedAsset.asset_digest)
            .join(AssetGroup, AssetGroup.id == Asset.asset_group_id)
            .where(AssetGroup.genome_digest == genome_digest)
        )
    return list(session.exec(statement.order_by(*_ORDER)).all())


def list_page(
    session: Session,
    *,
    asset_digest: str | None = None,
    mode: str | None = None,
    search: Search | None = None,
    offset: int,
    limit: int,
) -> tuple[list[StagedAsset], int]:
    """
    One page of staged rows, by asset digest and mode, and the total.

    Raises:
        ValueError: If the search names a field not in :data:`SEARCH_FIELDS`.
    """
    statement = select(StagedAsset)
    if asset_digest:
        statement = statement.where(StagedAsset.asset_digest == asset_digest)
    if mode:
        statement = statement.where(StagedAsset.mode == mode)
    if search is not None:
        condition = search_condition(
            StagedAsset, search.term, search.operator, search.fields, SEARCH_FIELDS
        )
        if condition is not None:
            statement = statement.where(condition)
    return paginate(session, statement, order_by=_ORDER, offset=offset, limit=limit)
