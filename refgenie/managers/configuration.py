from collections.abc import Iterable


from refgenie.db.queries import configuration as queries
from refgenie.db.tables import Configuration
from refgenie.managers.base import ResourceManager
from refgenie.models import Page


class ConfigurationManager(ResourceManager):
    """
    The configuration row: genome folders, servers, and the other settings kept
    in the database. Push remotes are ``RemoteManager`` (``rgc.remote``).
    """

    def list_all(self) -> Iterable[Configuration]:
        """
        List the configurations.

        Returns:
            list[Configuration]: The list of configurations.
        """
        with self._database_session as session:
            return queries.list_all(session)

    def get(self, id: int) -> Configuration:
        """
        Get a configuration by its ID.

        Args:
            id: The ID of the configuration.

        Returns:
            Configuration: The configuration.
        """
        with self._database_session as session:
            return queries.get(session, id)

    def get_latest(self) -> Configuration:
        """
        Get the latest configuration.

        Returns:
            Configuration: The latest configuration.
        """
        with self._database_session as session:
            return queries.get_latest(session)

    def page(self, *, offset: int, limit: int) -> Page[Configuration]:
        """One page of configuration rows, by id."""
        with self._database_session as session:
            rows, total = queries.list_page(session, offset=offset, limit=limit)
        return Page(items=rows, offset=offset, limit=limit, total=total)
