"""Queries on the ``configuration`` table."""

from sqlmodel import Session, desc, select

from refgenie.db.queries.common import one_or_raise
from refgenie.db.queries.paging import paginate
from refgenie.db.tables import Configuration
from refgenie.exceptions import MissingConfigDataError


def get_latest(session: Session) -> Configuration:
    """
    The configuration row in force: the highest id wins.

    Never assume id == 1; a restored dump or a merged database can leave the
    only row with any id.

    Raises:
        MissingConfigDataError: If there is no configuration row.
    """
    return one_or_raise(
        session,
        select(Configuration).order_by(desc(Configuration.id)).limit(1),
        MissingConfigDataError("No configuration row in the database"),
    )


def exists_any(session: Session) -> bool:
    """Whether any configuration row exists."""
    return session.exec(select(Configuration.id)).first() is not None


def get(session: Session, id: int) -> Configuration:
    """
    The configuration with this id.

    Raises:
        MissingConfigDataError: If there is none.
    """
    return one_or_raise(
        session,
        select(Configuration).where(Configuration.id == id),
        MissingConfigDataError(f"No configuration with id {id}"),
    )


def list_all(session: Session) -> list[Configuration]:
    """Every configuration row, by id."""
    return list(session.exec(select(Configuration).order_by(Configuration.id)).all())


def list_page(session: Session, *, offset: int, limit: int) -> tuple[list[Configuration], int]:
    """One page of configuration rows, by id, and the total."""
    return paginate(
        session, select(Configuration), order_by=(Configuration.id,), offset=offset, limit=limit
    )
