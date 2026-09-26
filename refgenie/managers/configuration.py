from collections.abc import Iterable

from rich.table import Table
from sqlmodel import select, desc

from refgenie.config import config
from refgenie.db.tables import Configuration
from refgenie.exceptions import MissingConfigDataError
from refgenie.managers.base import ResourceManager
from refgenie.managers.queries import one_or_raise
from refgenie.utils.tables import build_table


def latest_configuration(session) -> Configuration:
    """
    The configuration row in force, read through a caller's session.

    Same rule as `ConfigurationManager.get_latest`: the highest id wins. Never
    assume id == 1; a restored dump or a merged database can leave the only row
    with any id.
    """
    return one_or_raise(
        session,
        select(Configuration).order_by(desc(Configuration.id)),
        MissingConfigDataError("No configuration row in the database"),
        unique=True,
    )


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

        with self._database_session as session:
            configuration = session.exec(select(Configuration)).unique().one_or_none()
            if configuration is None:
                return None
        configuration = configuration.model_dump()
        configuration.pop("id")
        configuration.pop("updated_at")
        configuration.pop("created_at")

        # Format configuration values, special handling for servers list
        formatted_config = {}
        for k, v in configuration.items():
            if k == "servers" and isinstance(v, list):
                # Format servers as bulleted list
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
            return session.exec(select(Configuration)).unique().all()

    def get(self, id: int) -> Configuration:
        """
        Get a configuration by its ID.

        Args:
            id: The ID of the configuration.

        Returns:
            Configuration: The configuration.
        """
        with self._database_session as session:
            return one_or_raise(
                session,
                select(Configuration).where(Configuration.id == id),
                MissingConfigDataError(f"No configuration with id {id}"),
                unique=True,
            )

    def get_latest(self) -> Configuration:
        """
        Get the latest configuration.

        Returns:
            Configuration: The latest configuration.
        """
        with self._database_session as session:
            return one_or_raise(
                session,
                select(Configuration).order_by(desc(Configuration.id)),
                MissingConfigDataError("No configuration row in the database"),
                unique=True,
            )
