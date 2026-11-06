"""Queries on the ``store`` federation registry."""

from sqlmodel import Session, select

from refgenie.db.queries.common import one_or_raise
from refgenie.db.tables import Store
from refgenie.exceptions import MissingStoreError


def get(session: Session, name: str) -> Store:
    """
    The store row named ``name``.

    Raises:
        MissingStoreError: If there is none.
    """
    return one_or_raise(session, select(Store).where(Store.name == name), MissingStoreError(name))


def find(session: Session, name: str) -> Store | None:
    """The store row named ``name``, or None."""
    return session.exec(select(Store).where(Store.name == name)).one_or_none()


def list_all(session: Session, *, enabled_only: bool = False) -> list[Store]:
    """Store rows in priority order (lowest integer first), then by name."""
    statement = select(Store)
    if enabled_only:
        statement = statement.where(Store.enabled)
    return list(session.exec(statement.order_by(Store.priority, Store.name)).all())
