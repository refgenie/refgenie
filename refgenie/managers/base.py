from contextlib import contextmanager
from collections.abc import Generator

from sqlalchemy.engine import Engine
from sqlmodel import Session


class ResourceManager:
    """
    A base class for the resource manager to
    be used by Refgenie to manage independent databse resources.
    """

    def __init__(self, database_engine: Engine):
        self.database_engine = database_engine

    @property
    @contextmanager
    def _database_session(self) -> Generator[Session, None, None]:
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
        ``AssetManager._exists_in_session``) rather than calling a public
        manager method that opens its own.
        """
        session = Session(self.database_engine, expire_on_commit=False)
        try:
            yield session
        except:
            session.rollback()
            raise
        finally:
            session.close()
