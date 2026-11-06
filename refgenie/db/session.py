"""The session lifecycle every manager's ``_database_session`` uses."""

from collections.abc import Generator
from contextlib import contextmanager

from sqlalchemy.engine import Engine
from sqlalchemy.exc import IntegrityError
from sqlmodel import Session

from refgenie.exceptions import CatalogConflictError

__all__ = ["Session", "commit_or_conflict", "database_session"]


@contextmanager
def database_session(database_engine: Engine) -> Generator[Session, None, None]:
    """
    A new session on ``database_engine``, rolled back on error and always closed.

    ``expire_on_commit=False`` keeps rows readable after the session closes, so
    a db function's eager loads are all a caller needs.
    """
    session = Session(database_engine, expire_on_commit=False)
    try:
        yield session
    except:
        session.rollback()
        raise
    finally:
        session.close()


def commit_or_conflict(session: Session) -> None:
    """
    Commit, turning a uniqueness or key clash into ``CatalogConflictError``.

    Refgenie does no locking, so racing writers meet at the constraint. The
    session is rolled back first, so the caller can read the winner's row.
    """
    try:
        session.commit()
    except IntegrityError as e:
        session.rollback()
        raise CatalogConflictError(str(e.orig)) from e
