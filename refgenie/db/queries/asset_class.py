"""Queries on asset classes. Every asset class read loads its seek keys."""

from sqlalchemy.orm import selectinload
from sqlalchemy.sql.expression import null
from sqlmodel import Session, select

from refgenie.db.queries.common import one_or_raise, resolve_latest
from refgenie.db.queries.paging import paginate
from refgenie.db.queries.search import search_condition
from refgenie.db.tables import AssetClass, Recipe
from refgenie.exceptions import MissingAssetClassError
from refgenie.utils.search import Search

#: Fields an asset class search may match.
SEARCH_FIELDS = ["name", "version", "description"]

_ORDER = (AssetClass.name, AssetClass.version, AssetClass.id)


def _select(*, with_groups: bool = False):
    statement = select(AssetClass).options(selectinload(AssetClass.seek_keys))
    if with_groups:
        statement = statement.options(selectinload(AssetClass.asset_groups))
    return statement


def get(
    session: Session, name: str, version: str | None = None, *, with_groups: bool = False
) -> AssetClass:
    """
    The asset class ``name`` at ``version``, or its latest version by semver.

    Args:
        session: The session to execute in.
        name: The asset class name.
        version: The version, or None for the latest.
        with_groups: Also load the asset groups of this class.

    Raises:
        MissingAssetClassError: If there is no such asset class.
    """
    label = f"{name} v{version}" if version else name
    if version is None:
        version = resolve_latest(session, AssetClass, name, MissingAssetClassError(label)).version
    return one_or_raise(
        session,
        _select(with_groups=with_groups).where(
            AssetClass.name == name, AssetClass.version == version
        ),
        MissingAssetClassError(label),
        unique=True,
    )


def get_by_id(session: Session, id: int) -> AssetClass | None:
    """The asset class with this id, or None."""
    return session.exec(_select().where(AssetClass.id == id)).unique().one_or_none()


def list_versions(session: Session, name: str, *, with_groups: bool = False) -> list[AssetClass]:
    """Every registered version of the asset class ``name``."""
    statement = _select(with_groups=with_groups).where(AssetClass.name == name).order_by(*_ORDER)
    return list(session.exec(statement).unique().all())


def exists(session: Session, name: str, version: str) -> bool:
    """Whether ``name`` at ``version`` is registered."""
    statement = select(AssetClass.id).where(AssetClass.name == name, AssetClass.version == version)
    return session.exec(statement).first() is not None


def exists_any_version(session: Session, name: str) -> bool:
    """Whether any version of ``name`` is registered."""
    return session.exec(select(AssetClass.id).where(AssetClass.name == name)).first() is not None


def list_all(session: Session) -> list[AssetClass]:
    """Every asset class, by name and version."""
    return list(session.exec(_select().order_by(*_ORDER)).unique().all())


def list_root(session: Session) -> list[AssetClass]:
    """Asset classes some recipe builds with no input asset classes."""
    statement = _select().join(Recipe).where(Recipe.input_asset_classes == null()).order_by(*_ORDER)
    return list(session.exec(statement).unique().all())


def list_page(
    session: Session,
    *,
    name: str | None = None,
    version: str | None = None,
    search: Search | None = None,
    offset: int,
    limit: int,
) -> tuple[list[AssetClass], int]:
    """
    One page of asset classes, by name and version, and the total.

    Raises:
        ValueError: If the search names a field not in :data:`SEARCH_FIELDS`.
    """
    statement = _select()
    if name:
        statement = statement.where(AssetClass.name == name)
    if version:
        statement = statement.where(AssetClass.version == version)
    if search is not None:
        condition = search_condition(
            AssetClass, search.term, search.operator, search.fields, SEARCH_FIELDS
        )
        if condition is not None:
            statement = statement.where(condition)
    return paginate(session, statement, order_by=_ORDER, offset=offset, limit=limit)
