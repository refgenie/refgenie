"""
Executing a lookup that must find exactly one row.

"Look this up, and raise a refgenie exception if it isn't there" is the single
most common shape in the managers, and it used to be hand-rolled three
different ways: catching :class:`NoResultFound` around ``.one()``, testing
``.one_or_none()`` against ``None``, or letting the raw SQLAlchemy exception
escape into a caller that only knows about `RefgenieError`. `one_or_raise`
gives all of them one spelling and one failure mode.

Keep this module dependent only on ``refgenie.db`` and SQLModel/SQLAlchemy --
see ``refgenie.managers.asset.queries``, which holds the SELECT builders this
executes, for the same rule.
"""

from typing import TypeVar

from sqlalchemy.exc import NoResultFound
from sqlmodel import Session
from sqlmodel.sql.expression import SelectOfScalar

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
