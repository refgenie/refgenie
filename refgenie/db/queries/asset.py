"""
Queries on assets, their names, their seek keys and their parent/child links.

Every function that returns an ``Asset`` loads the same set: its group with the
group's genome and asset class (``serving_modes``, ``asset_class_name``,
``registry_path``), its name rows and its seek keys. ``with_links=True`` also
loads its parents and children, each with that same set.
"""

from collections.abc import Iterable

from sqlalchemy.orm import selectinload
from sqlmodel import Session, func, select

from refgenie.db.queries.common import one_or_raise
from refgenie.db.queries.paging import paginate
from refgenie.db.queries.search import search_condition
from refgenie.db.tables import Asset, AssetGroup, AssetName, Recipe, SeekKey
from refgenie.exceptions import MissingAssetError, MissingSeekKeyError
from refgenie.utils.search import Search

#: Fields an asset search may match.
SEARCH_FIELDS = ["name", "digest", "path"]

#: By genome, then group name, then digest (unique).
_ORDER = (AssetGroup.genome_digest, AssetGroup.name, Asset.digest)


def _loads(relationship=None) -> list:
    """The loader options for an asset (or, given a relationship, its relatives)."""
    if relationship is None:
        return [
            selectinload(Asset.asset_group).selectinload(AssetGroup.genome),
            selectinload(Asset.asset_group).selectinload(AssetGroup.asset_class),
            selectinload(Asset.asset_names),
            selectinload(Asset.seek_keys),
        ]
    return [
        selectinload(relationship).selectinload(Asset.asset_group).selectinload(AssetGroup.genome),
        selectinload(relationship)
        .selectinload(Asset.asset_group)
        .selectinload(AssetGroup.asset_class),
        selectinload(relationship).selectinload(Asset.asset_names),
        selectinload(relationship).selectinload(Asset.seek_keys),
    ]


def _select(*, with_links: bool = False):
    """``select(Asset)`` with the standard loads, and the links' if asked."""
    options = _loads()
    if with_links:
        options += _loads(Asset.parents) + _loads(Asset.children)
    return select(Asset).options(*options)


def _by_name(statement, genome_digest: str, asset_group_name: str, asset_name: str):
    """Restrict ``statement`` (over Asset) to ``genome_digest/asset_group_name:asset_name``."""
    return (
        statement.join(AssetName, AssetName.asset_digest == Asset.digest)
        .join(AssetGroup, AssetGroup.id == AssetName.asset_group_id)
        .where(
            AssetName.name == asset_name,
            AssetGroup.name == asset_group_name,
            AssetGroup.genome_digest == genome_digest,
        )
    )


def find_by_digest(session: Session, digest: str, *, with_links: bool = False) -> Asset | None:
    """The asset holding content ``digest``, or None."""
    statement = _select(with_links=with_links).where(Asset.digest == digest)
    return session.exec(statement).unique().one_or_none()


def get_by_digest(session: Session, digest: str, *, with_links: bool = False) -> Asset:
    """
    The asset holding content ``digest``.

    Raises:
        MissingAssetError: If there is none.
    """
    return one_or_raise(
        session,
        _select(with_links=with_links).where(Asset.digest == digest),
        MissingAssetError(digest=digest),
        unique=True,
    )


def find_by_name(
    session: Session, genome_digest: str, asset_group_name: str, asset_name: str
) -> Asset | None:
    """The asset named ``genome_digest/asset_group_name:asset_name``, or None."""
    statement = _by_name(_select(), genome_digest, asset_group_name, asset_name)
    return session.exec(statement).unique().one_or_none()


def get_by_name(
    session: Session,
    genome_digest: str,
    asset_group_name: str,
    asset_name: str,
    *,
    with_links: bool = False,
) -> Asset:
    """
    The asset named ``genome_digest/asset_group_name:asset_name``.

    Raises:
        MissingAssetError: If no asset has that name.
    """
    return one_or_raise(
        session,
        _by_name(_select(with_links=with_links), genome_digest, asset_group_name, asset_name),
        MissingAssetError(genome=genome_digest, asset_group=asset_group_name, asset=asset_name),
        unique=True,
    )


def exists_name(
    session: Session, genome_digest: str, asset_group_name: str, asset_name: str
) -> bool:
    """Whether a name is taken in a group. False for an unknown genome."""
    return find_name(session, genome_digest, asset_group_name, asset_name) is not None


def find_name(
    session: Session, genome_digest: str, asset_group_name: str, asset_name: str
) -> AssetName | None:
    """The ``AssetName`` row for ``genome_digest/asset_group_name:asset_name``, or None."""
    return session.exec(
        select(AssetName)
        .join(AssetGroup, AssetGroup.id == AssetName.asset_group_id)
        .where(
            AssetName.name == asset_name,
            AssetGroup.name == asset_group_name,
            AssetGroup.genome_digest == genome_digest,
        )
    ).one_or_none()


def list_all(
    session: Session,
    genome_digests: Iterable[str] | None = None,
    asset_group_name: str | None = None,
) -> list[Asset]:
    """
    Assets, by genome, group and digest.

    Args:
        session: The session to execute in.
        genome_digests: Only these genomes. None reads every asset.
        asset_group_name: Only groups with this name.
    """
    statement = _select().join(AssetGroup, AssetGroup.id == Asset.asset_group_id)
    if genome_digests is not None:
        statement = statement.where(AssetGroup.genome_digest.in_(list(genome_digests)))
    if asset_group_name is not None:
        statement = statement.where(AssetGroup.name == asset_group_name)
    return list(session.exec(statement.order_by(*_ORDER)).unique().all())


def list_by_digests(session: Session, digests: Iterable[str]) -> list[Asset]:
    """The assets holding these content digests, by digest. Unknown ones are skipped."""
    statement = _select().where(Asset.digest.in_(list(digests))).order_by(Asset.digest)
    return list(session.exec(statement).unique().all())


def list_page(
    session: Session,
    *,
    genome_digest: str | None = None,
    asset_group_name: str | None = None,
    asset_group_id: int | None = None,
    name: str | None = None,
    recipe_name: str | None = None,
    search: Search | None = None,
    offset: int,
    limit: int,
) -> tuple[list[Asset], int]:
    """
    One page of assets, by genome, group and digest, and the total.

    ``name`` matches any of an asset's names, not only the canonical one,
    which is what lets a client pull by a non-canonical name.

    Raises:
        ValueError: If the search names a field not in :data:`SEARCH_FIELDS`.
    """
    # Join assetgroup once, with an explicit ON: assetname gives a bare join
    # two foreign-key paths.
    statement = _select().join(AssetGroup, AssetGroup.id == Asset.asset_group_id)
    if name:
        statement = statement.join(AssetName, AssetName.asset_digest == Asset.digest).where(
            AssetName.name == name
        )
    if asset_group_id is not None:
        statement = statement.where(AssetGroup.id == asset_group_id)
    if asset_group_name:
        statement = statement.where(AssetGroup.name == asset_group_name)
    if genome_digest:
        statement = statement.where(AssetGroup.genome_digest == genome_digest)
    if recipe_name:
        statement = statement.join(Recipe, Recipe.id == Asset.recipe_id).where(
            Recipe.name == recipe_name
        )
    if search is not None:
        condition = search_condition(
            Asset, search.term, search.operator, search.fields, SEARCH_FIELDS
        )
        if condition is not None:
            statement = statement.where(condition)
    return paginate(session, statement, order_by=_ORDER, offset=offset, limit=limit)


def get_seek_key(
    session: Session,
    genome_digest: str,
    asset_group_name: str,
    asset_name: str,
    seek_key_name: str,
) -> SeekKey:
    """
    One named seek key of ``genome_digest/asset_group_name:asset_name``, its asset loaded.

    Raises:
        MissingSeekKeyError: If there is none.
    """
    statement = (
        select(SeekKey)
        .join(Asset)
        .join(AssetName, AssetName.asset_digest == Asset.digest)
        .join(AssetGroup, AssetGroup.id == AssetName.asset_group_id)
        .options(selectinload(SeekKey.asset))
        .where(
            SeekKey.name == seek_key_name,
            AssetName.name == asset_name,
            AssetGroup.name == asset_group_name,
            AssetGroup.genome_digest == genome_digest,
        )
    )
    return one_or_raise(
        session,
        statement,
        MissingSeekKeyError(
            genome=genome_digest,
            asset_group=asset_group_name,
            asset=asset_name,
            seek_key=seek_key_name,
        ),
        unique=True,
    )


def count(session: Session) -> int:
    """How many assets there are."""
    return session.exec(select(func.count()).select_from(Asset)).one()
