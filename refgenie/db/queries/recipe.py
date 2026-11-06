"""
Queries on recipes. Every recipe read loads its output asset class (with seek
keys) and its input asset class links (with their asset classes).
"""

from sqlalchemy.orm import selectinload
from sqlmodel import Session, select

from refgenie.db.queries.common import one_or_raise, resolve_latest
from refgenie.db.queries.paging import paginate
from refgenie.db.queries.search import search_condition
from refgenie.db.tables import Asset, AssetClass, Recipe, RecipeAssetClassesInputs
from refgenie.exceptions import MissingRecipeError
from refgenie.utils.search import Search

#: Fields a recipe search may match.
SEARCH_FIELDS = ["name", "version", "description"]

_ORDER = (Recipe.name, Recipe.version, Recipe.id)


def _select():
    return select(Recipe).options(
        selectinload(Recipe.output_asset_class).selectinload(AssetClass.seek_keys),
        selectinload(Recipe.input_asset_classes).selectinload(RecipeAssetClassesInputs.asset_class),
    )


def get(session: Session, name: str, version: str | None = None) -> Recipe:
    """
    The recipe ``name`` at ``version``, or its latest version by semver.

    Raises:
        MissingRecipeError: If there is no such recipe.
    """
    label = f"{name} v{version}" if version else name
    if version is None:
        version = resolve_latest(session, Recipe, name, MissingRecipeError(label)).version
    return one_or_raise(
        session,
        _select().where(Recipe.name == name, Recipe.version == version),
        MissingRecipeError(label),
        unique=True,
    )


def get_by_id(session: Session, id: int) -> Recipe | None:
    """The recipe with this id, or None."""
    return session.exec(_select().where(Recipe.id == id)).unique().one_or_none()


def exists(session: Session, name: str, version: str | None = None) -> bool:
    """Whether recipe ``name`` (at ``version``, if given) is registered."""
    statement = select(Recipe.id).where(Recipe.name == name)
    if version:
        statement = statement.where(Recipe.version == version)
    return session.exec(statement).first() is not None


def list_versions(session: Session, name: str) -> list[Recipe]:
    """Every registered version of recipe ``name``."""
    return list(session.exec(_select().where(Recipe.name == name).order_by(*_ORDER)).unique().all())


def list_all(session: Session) -> list[Recipe]:
    """Every recipe, by name and version."""
    return list(session.exec(_select().order_by(*_ORDER)).unique().all())


def list_by_output_asset_class(session: Session, asset_class_name: str) -> list[Recipe]:
    """Recipes producing any version of the asset class ``asset_class_name``."""
    statement = (
        _select()
        .join(AssetClass, AssetClass.id == Recipe.output_asset_class_id)
        .where(AssetClass.name == asset_class_name)
        .order_by(*_ORDER)
    )
    return list(session.exec(statement).unique().all())


def find_using_asset_class(session: Session, asset_class_id: int) -> Recipe | None:
    """A recipe that outputs or takes as input this asset class, or None."""
    output = session.exec(
        select(Recipe).where(Recipe.output_asset_class_id == asset_class_id).order_by(Recipe.id)
    ).first()
    if output is not None:
        return output
    return session.exec(
        select(Recipe)
        .join(RecipeAssetClassesInputs)
        .where(RecipeAssetClassesInputs.asset_class_id == asset_class_id)
        .order_by(Recipe.id)
    ).first()


def list_input_asset_classes(session: Session, recipe_id: int) -> list[tuple[AssetClass, str]]:
    """``(asset class, default)`` for each input asset of a recipe."""
    statement = (
        select(AssetClass, RecipeAssetClassesInputs.default)
        .join(RecipeAssetClassesInputs, RecipeAssetClassesInputs.asset_class_id == AssetClass.id)
        .where(RecipeAssetClassesInputs.recipe_id == recipe_id)
        .order_by(RecipeAssetClassesInputs.id)
    )
    return list(session.exec(statement).unique().all())


def list_assets_built(session: Session, recipe_id: int) -> list[Asset]:
    """The assets built with this recipe, by digest."""
    statement = select(Asset).where(Asset.recipe_id == recipe_id).order_by(Asset.digest)
    return list(session.exec(statement).unique().all())


def list_page(
    session: Session,
    *,
    output_asset_class: str | None = None,
    name: str | None = None,
    version: str | None = None,
    search: Search | None = None,
    offset: int,
    limit: int,
) -> tuple[list[Recipe], int]:
    """
    One page of recipes, by name and version, and the total.

    Raises:
        ValueError: If the search names a field not in :data:`SEARCH_FIELDS`.
    """
    statement = _select()
    if output_asset_class:
        statement = statement.join(AssetClass, AssetClass.id == Recipe.output_asset_class_id).where(
            AssetClass.name == output_asset_class
        )
    if name:
        statement = statement.where(Recipe.name == name)
    if version:
        statement = statement.where(Recipe.version == version)
    if search is not None:
        condition = search_condition(
            Recipe, search.term, search.operator, search.fields, SEARCH_FIELDS
        )
        if condition is not None:
            statement = statement.where(condition)
    return paginate(session, statement, order_by=_ORDER, offset=offset, limit=limit)
