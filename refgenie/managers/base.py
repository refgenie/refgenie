from contextlib import AbstractContextManager, contextmanager
from collections.abc import Generator

from sqlalchemy.engine import Engine
from sqlmodel import Session


@contextmanager
def database_session(database_engine: Engine) -> Generator[Session, None, None]:
    """
    A new session on ``database_engine``, rolled back on error and always closed.

    Shared by every ``_database_session`` property: each manager's, and the
    root ``Refgenie``'s.
    """
    session = Session(database_engine, expire_on_commit=False)
    try:
        yield session
    except:
        session.rollback()
        raise
    finally:
        session.close()


class ResourceManager:
    """
    A base class for the resource manager to
    be used by Refgenie to manage independent databse resources.
    """

    def __init__(self, database_engine: Engine):
        self.database_engine = database_engine

    @property
    def _database_session(self) -> AbstractContextManager[Session]:
        """
        Provide a transactional scope around a series of query
        operations.

        This is a property, so *every access constructs a new Session* against
        the same engine. It must therefore not be re-entered from inside an open
        block: doing so opens a second connection while the first still holds an
        uncommitted transaction. Under SQLite that is a write-lock hazard, and
        under any backend the inner session reads a snapshot that cannot see the
        outer one's pending writes. A method that needs to query while a session
        is open must take that session as an argument (see
        ``refgenie.managers.asset.content.exists_in_session``) rather than calling a public
        manager method that opens its own.
        """
        return database_session(self.database_engine)
