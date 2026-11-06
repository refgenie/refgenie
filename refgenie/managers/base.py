from contextlib import AbstractContextManager

from refgenie.db.engine import Engine
from refgenie.db.session import Session, database_session


class ResourceManager:
    """
    Base class for managers of database resources.
    """

    def __init__(self, database_engine: Engine):
        self.database_engine = database_engine

    @property
    def _database_session(self) -> AbstractContextManager[Session]:
        """
        Provide a transactional scope around a series of query operations.

        Every access opens a new Session, so do not nest it: an inner session
        cannot see the outer one's pending writes, and under SQLite it risks a
        write-lock. Inside an open block, pass the session to
        ``refgenie.db.queries`` functions instead.
        """
        return database_session(self.database_engine)
