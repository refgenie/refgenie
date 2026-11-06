"""
Lookups that must find exactly one row, raising a caller-built refgenie
exception instead of SQLAlchemy's :class:`NoResultFound`.
"""

from typing import TypeVar

from sqlalchemy.exc import NoResultFound
from sqlmodel import Session, select
from sqlmodel.sql.expression import SelectOfScalar

from refgenie.utils.versioning import semver_sort_key

T = TypeVar("T")


def one_or_raise(
    session: Session,
    statement: SelectOfScalar[T],
    exc: Exception,
    *,
    unique: bool = False,
) -> T:
    """
    Execute ``statement`` and return its single row, or raise ``exc``.

    Args:
        session: The session to execute in.
        statement: A statement selecting exactly one row.
        exc: The exception to raise when nothing matched.
        unique: Apply ``.unique()`` before reading the row. Required when the
            statement joined-eager-loads a collection, which otherwise yields
            the parent row once per child.

    Returns:
        The single matching row.

    Raises:
        Exception: ``exc``, if nothing matched.
        MultipleResultsFound: If more than one distinct row matched.
    """
    result = session.exec(statement)
    if unique:
        result = result.unique()
    try:
        return result.one()
    except NoResultFound:
        raise exc from None


def resolve_latest(session: Session, model_class, name: str, exc: Exception):
    """Query the DB for the latest version of a given name using semver ordering.

    Works for any SQLModel with `name` and `version` fields (AssetClass, Recipe).

    Args:
        session: The database session.
        model_class: The SQLModel class (e.g., AssetClass or Recipe).
        name: The name to look up.
        exc: The exception to raise when no row has that name.

    Returns:
        The model instance with the latest version.

    Raises:
        Exception: ``exc``, if no row has that name.
    """
    results = session.exec(select(model_class).where(model_class.name == name)).unique().all()

    if not results:
        raise exc

    if len(results) == 1:
        return results[0]

    return max(results, key=lambda r: semver_sort_key(r.version))
