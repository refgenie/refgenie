"""
Executing a lookup that must find exactly one row.

"Look this up, and raise a refgenie exception if it isn't there" is the single
most common shape in the managers. `one_or_raise` is its one spelling: it
returns the single matching row, or raises the exception the caller built
(normally a `RefgenieError` naming what was looked up), so a raw SQLAlchemy
:class:`NoResultFound` never escapes into a caller that only knows about
`RefgenieError`. Do not hand-roll this with ``.one()`` or
``.one_or_none()`` in a manager. `resolve_latest` is the
versioned variant: the highest-semver row for a name.

Keep this module dependent only on ``refgenie.db`` and SQLModel/SQLAlchemy --
see ``refgenie.managers.asset.queries``, which holds the SELECT builders this
executes, for the same rule.
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
        exc: The exception to raise when nothing matched. Built by the caller,
            so it carries the names the caller was looking for.
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


def resolve_latest(session: Session, model_class, name: str):
    """Query the DB for the latest version of a given name using semver ordering.

    Works for any SQLModel with `name` and `version` fields (AssetClass, Recipe).

    Args:
        session: The database session.
        model_class: The SQLModel class (e.g., AssetClass or Recipe).
        name: The name to look up.

    Returns:
        The model instance with the latest version.

    Raises:
        NoResultFound: If no entries exist with that name.
    """
    results = session.exec(select(model_class).where(model_class.name == name)).unique().all()

    if not results:
        raise NoResultFound(f"No {model_class.__name__} found with name '{name}'")

    if len(results) == 1:
        return results[0]

    return max(results, key=lambda r: semver_sort_key(r.version))
