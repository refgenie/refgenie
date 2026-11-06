from collections.abc import Iterable

from rich.table import Table

from refgenie.config import config
from refgenie.db.queries import configuration as queries
from refgenie.db.tables import Configuration
from refgenie.exceptions import MissingConfigDataError
from refgenie.managers.base import ResourceManager
from refgenie.models import Page
from refgenie.utils.tables import build_table


class ConfigurationManager(ResourceManager):
    """
    The configuration row: genome folders, servers, and the other settings kept
    in the database. Push remotes are ``RemoteManager`` (``rgc.remote``).
    """

    def table(self) -> Table | None:
        """
        Return the configuration table.

        Returns:
            Table: The configuration table.
        """

        try:
            configuration = self.get_latest()
        except MissingConfigDataError:
            return None
        configuration = configuration.model_dump()
        configuration.pop("id")
        configuration.pop("updated_at")
        configuration.pop("created_at")

        formatted_config = {}
        for k, v in configuration.items():
            if k == "servers" and isinstance(v, list):
                formatted_config[k] = "\n".join(f"• {server}" for server in v)
            else:
                formatted_config[k] = str(v)
        configuration = formatted_config
        return build_table(
            "Refgenie configuration",
            list(configuration.keys()),
            [list(configuration.values())],
            caption="\n".join(
                ["\nEnvironment-based configuration\n"]
                + [f"{k}: {v}" for k, v in config.model_dump().items()]
            ),
            caption_justify="left",
            caption_style="italic",
        )

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
