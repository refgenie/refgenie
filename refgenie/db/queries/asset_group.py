"""Queries on asset groups and the default-name flag of their asset names."""

from sqlalchemy import update
from sqlalchemy.orm import selectinload
from sqlmodel import Session, select

from refgenie.db.queries.common import one_or_raise
from refgenie.db.queries.paging import paginate
from refgenie.db.queries.search import search_condition
from refgenie.db.tables import AssetClass, AssetGroup, AssetName
from refgenie.exceptions import MissingAssetGroupError
from refgenie.utils.search import Search

#: Fields an asset group search may match.
SEARCH_FIELDS = ["name"]

_ORDER = (AssetGroup.genome_digest, AssetGroup.name, AssetGroup.id)


def get(session: Session, genome_digest: str, asset_group_name: str) -> AssetGroup:
    """
    The asset group named ``genome_digest/asset_group_name``.

    Raises:
        MissingAssetGroupError: If there is none.
    """
    return one_or_raise(
        session,
        select(AssetGroup).where(
            AssetGroup.name == asset_group_name, AssetGroup.genome_digest == genome_digest
        ),
        MissingAssetGroupError(genome=genome_digest, asset_group=asset_group_name),
    )


def get_by_id(session: Session, id: int) -> AssetGroup | None:
    """The asset group with this id, or None."""
    return session.get(AssetGroup, id)


def exists(session: Session, genome_digest: str, asset_group_name: str) -> bool:
    """Whether ``genome_digest/asset_group_name`` exists."""
    statement = select(AssetGroup.id).where(
        AssetGroup.name == asset_group_name, AssetGroup.genome_digest == genome_digest
    )
    return session.exec(statement).first() is not None


def exists_for_genome(session: Session, genome_digest: str) -> bool:
    """Whether the genome has any asset group left."""
    statement = select(AssetGroup.id).where(AssetGroup.genome_digest == genome_digest)
    return session.exec(statement).first() is not None


def list_all(session: Session, genome_digests: list[str] | None = None) -> list[AssetGroup]:
    """Asset groups with their genome and assets loaded, by genome then name."""
    statement = select(AssetGroup).options(
        selectinload(AssetGroup.genome), selectinload(AssetGroup.assets)
    )
    if genome_digests is not None:
        statement = statement.where(AssetGroup.genome_digest.in_(genome_digests))
    return list(session.exec(statement.order_by(*_ORDER)).unique().all())


def list_page(
    session: Session,
    *,
    genome_digest: str | None = None,
    asset_class: str | None = None,
    name: str | None = None,
    id: int | None = None,
    search: Search | None = None,
    offset: int,
    limit: int,
) -> tuple[list[AssetGroup], int]:
    """
    One page of asset groups, by genome, name and id, and the total.

    Raises:
        ValueError: If the search names a field not in :data:`SEARCH_FIELDS`.
    """
    statement = select(AssetGroup)
    if genome_digest:
        statement = statement.where(AssetGroup.genome_digest == genome_digest)
    if asset_class:
        statement = statement.join(AssetClass).where(AssetClass.name == asset_class)
    if name:
        statement = statement.where(AssetGroup.name == name)
    if id is not None:
        statement = statement.where(AssetGroup.id == id)
    if search is not None:
        condition = search_condition(
            AssetGroup, search.term, search.operator, search.fields, SEARCH_FIELDS
        )
        if condition is not None:
            statement = statement.where(condition)
    return paginate(session, statement, order_by=_ORDER, offset=offset, limit=limit)


def find_name(session: Session, asset_group_id: int, asset_name: str) -> AssetName | None:
    """The ``AssetName`` row for this name in this group, or None."""
    return session.exec(
        select(AssetName).where(
            AssetName.asset_group_id == asset_group_id, AssetName.name == asset_name
        )
    ).one_or_none()


def get_default_name(session: Session, asset_group_id: int) -> str | None:
    """The group's default asset name, or None if no default is set."""
    return session.exec(
        select(AssetName.name).where(
            AssetName.asset_group_id == asset_group_id,
            AssetName.is_default == True,  # noqa: E712 - SQL boolean, not Python identity
        )
    ).one_or_none()


def clear_default(session: Session, asset_group_id: int) -> None:
    """Unflag the group's current default name, if any."""
    session.exec(
        update(AssetName)
        .where(
            AssetName.asset_group_id == asset_group_id,
            AssetName.is_default == True,  # noqa: E712 - SQL boolean, not Python identity
        )
        .values(is_default=False)
    )
