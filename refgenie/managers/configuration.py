from collections.abc import Iterable

from rich.table import Table
from sqlalchemy.orm import selectinload
from sqlmodel import select, desc

from refgenie.config import config
from collections import defaultdict

from refgenie.db.tables import Configuration, Remote, RemoteAssetLink, RemoteType, StagedAsset
from refgenie.exceptions import MissingConfigDataError, MissingRemoteError
from refgenie.logger import logger
from refgenie.managers.base import ResourceManager
from refgenie.managers.queries import one_or_raise
from refgenie.utils.tables import build_table


class ConfigurationManager(ResourceManager):
    """
    A manager for configuration.
    """

    @staticmethod
    def _latest(session) -> Configuration:
        """
        The configuration row this manager operates on, read through a caller's session.

        Same rule as `get_latest`: the highest id wins. Never assume id == 1; a
        restored dump or a merged database can leave the only row with any id.
        """
        return one_or_raise(
            session,
            select(Configuration).order_by(desc(Configuration.id)),
            MissingConfigDataError("No configuration row in the database"),
            unique=True,
        )

    def subscribe(self, server_urls: str | list[str], reset: bool = False):
        """
        Subscribe to a list of servers.

        Args:
            server_urls: Server URL or list of server URLs to subscribe to.
            reset: If True, overwrite the current list instead of adding to it.
        """
        if isinstance(server_urls, str):
            server_urls = [server_urls]
        with self._database_session as session:
            configuration = self._latest(session)
            configuration.servers = (
                list(server_urls) if reset else list(set(configuration.servers).union(server_urls))
            )
            session.add(configuration)
            session.commit()
        logger.info(f"Subscribed to servers: {' '.join(server_urls)}")

    def unsubscribe(self, server_urls: list[str]):
        """
        Unsubscribe from a list of servers.

        Args:
            server_urls: The list of server URLs to unsubscribe from.
        """
        if not isinstance(server_urls, list):
            raise TypeError(f"servers must be a list of strings, not {type(server_urls)}")
        with self._database_session as session:
            configuration = self._latest(session)
            configuration.servers = list(set(configuration.servers).difference(server_urls))
            session.add(configuration)
            session.commit()
        logger.info(f"Unsubscribed from servers: {server_urls}")

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

    def remote_exists(self, remote_type: str | RemoteType) -> bool:
        """
        Check if any remote of the specified type exists.

        Args:
            remote_type: The type of the remote to check.

        Returns:
            bool: True if at least one remote of this type exists, False otherwise.
        """
        with self._database_session as session:
            return (
                session.exec(select(Remote).where(Remote.type == RemoteType(remote_type))).first()
                is not None
            )

    def add_remote(
        self,
        remote_type: RemoteType,
        prefix: str,
        description: str,
        push_command: str | None = None,
    ) -> Remote:
        """
        Add a remote. Multiple remotes of the same type are allowed.

        Args:
            remote_type: The type of the remote.
            prefix: The prefix of the remote.
            description: The description of the remote.
            push_command: Optional shell command template for pushing assets.

        Returns:
            Remote: The added remote.
        """
        latest_configuration = self.get_latest()
        with self._database_session as session:
            remote = Remote(
                type=remote_type,
                prefix=prefix,
                description=description,
                push_command=push_command,
                configuration_id=latest_configuration.id,
            )
            session.add(remote)
            session.commit()
            session.refresh(remote)
            logger.info(f"Added {remote}")
            return remote

    def upsert_remote(
        self,
        name: str,
        type: RemoteType,
        prefix: str,
        push_command: str | None = None,
    ) -> Remote:
        """
        Upsert a remote by name (description field).

        If a remote with the given name (description) exists, update its fields.
        Otherwise create a new one.

        Args:
            name: The name (stored in description field) used for matching.
            type: The remote type.
            prefix: The remote prefix URL.
            push_command: Optional shell command template for pushing assets.

        Returns:
            Remote: The upserted remote.
        """
        with self._database_session as session:
            existing = session.exec(select(Remote).where(Remote.description == name)).first()
            if existing:
                existing.type = type
                existing.prefix = prefix
                existing.push_command = push_command
                session.add(existing)
                session.commit()
                session.refresh(existing)
                logger.info(f"Updated remote '{name}' (id={existing.id})")
                return existing
            else:
                # Queried through this session: get_latest opens its own, and
                # _database_session must not be re-entered from inside an open
                # block.
                latest_configuration = self._latest(session)
                remote = Remote(
                    type=type,
                    prefix=prefix,
                    description=name,
                    push_command=push_command,
                    configuration_id=latest_configuration.id,
                )
                session.add(remote)
                session.commit()
                session.refresh(remote)
                logger.info(f"Created remote '{name}' (id={remote.id})")
                return remote

    def remove_remote(
        self, remote_id: int | None = None, remote_type: str | RemoteType | None = None
    ) -> None:
        """
        Remove a remote by its ID or type.

        Args:
            remote_id: The ID of the remote to remove.
            remote_type: The type of the remote to remove (removes first match).

        Raises:
            MissingRemoteError: If the remote is not found.
            ValueError: If neither remote_id nor remote_type is provided.
        """
        with self._database_session as session:
            if remote_id is not None:
                remote = session.exec(select(Remote).where(Remote.id == remote_id)).first()
            elif remote_type is not None:
                remote_type = RemoteType(remote_type)
                remote = session.exec(select(Remote).where(Remote.type == remote_type)).first()
            else:
                raise ValueError("Must provide either remote_id or remote_type")

            if remote is None:
                raise MissingRemoteError(remote_id or remote_type)

            session.delete(remote)
            session.commit()
            logger.info(f"Removed {remote}")

    def link_asset_to_remote(
        self,
        remote_id: int,
        asset_digest: str,
        mode: str,
        pushed: bool = False,
    ) -> RemoteAssetLink:
        """
        Create a RemoteAssetLink record linking a staged asset to a remote.

        Args:
            remote_id: The ID of the remote.
            asset_digest: The asset digest.
            mode: The staging mode ("archive" or "file").
            pushed: Whether the asset has been pushed (default False = intent).

        Returns:
            RemoteAssetLink: The created link record.

        Raises:
            ValueError: If the (asset_digest, mode) StagedAsset does not exist.
        """
        with self._database_session as session:
            one_or_raise(
                session,
                select(StagedAsset).where(
                    StagedAsset.asset_digest == asset_digest,
                    StagedAsset.mode == mode,
                ),
                ValueError(
                    f"StagedAsset not found for {asset_digest=}, {mode=}. "
                    f"Asset must be staged before linking to a remote."
                ),
            )

            link = RemoteAssetLink(
                remote_id=remote_id,
                asset_digest=asset_digest,
                mode=mode,
                pushed=pushed,
            )
            session.add(link)
            session.commit()
            session.refresh(link)
            logger.info(f"Linked {asset_digest} ({mode}) to remote {remote_id} (pushed={pushed})")
            return link

    def unlink_asset_from_remote(
        self,
        remote_id: int,
        asset_digest: str,
        mode: str,
    ) -> None:
        """
        Remove a RemoteAssetLink record.

        Args:
            remote_id: The ID of the remote.
            asset_digest: The asset digest.
            mode: The staging mode.
        """
        with self._database_session as session:
            link = one_or_raise(
                session,
                select(RemoteAssetLink).where(
                    RemoteAssetLink.remote_id == remote_id,
                    RemoteAssetLink.asset_digest == asset_digest,
                    RemoteAssetLink.mode == mode,
                ),
                ValueError(
                    f"No link found for remote_id={remote_id}, "
                    f"asset_digest={asset_digest}, mode={mode}"
                ),
            )
            session.delete(link)
            session.commit()
            logger.info(f"Unlinked {asset_digest} ({mode}) from remote {remote_id}")

    def mark_pushed(
        self,
        remote_id: int,
        asset_digest: str,
        mode: str,
    ) -> None:
        """
        Set pushed=True on an existing RemoteAssetLink.

        Called by refgenie push after successful upload.

        Args:
            remote_id: The ID of the remote.
            asset_digest: The asset digest.
            mode: The staging mode.
        """
        with self._database_session as session:
            link = one_or_raise(
                session,
                select(RemoteAssetLink).where(
                    RemoteAssetLink.remote_id == remote_id,
                    RemoteAssetLink.asset_digest == asset_digest,
                    RemoteAssetLink.mode == mode,
                ),
                ValueError(
                    f"No link found for remote_id={remote_id}, "
                    f"asset_digest={asset_digest}, mode={mode}"
                ),
            )
            link.pushed = True
            session.add(link)
            session.commit()
            logger.info(f"Marked pushed: {asset_digest} ({mode}) on remote {remote_id}")

    def get_unpushed_links(
        self,
        remote_id: int | None = None,
    ) -> list[RemoteAssetLink]:
        """
        Query for RemoteAssetLink records where pushed=False.

        Args:
            remote_id: Optional filter by remote ID.

        Returns:
            List of unpushed RemoteAssetLink records.
        """
        with self._database_session as session:
            query = select(RemoteAssetLink).where(RemoteAssetLink.pushed == False)  # noqa: E712
            if remote_id is not None:
                query = query.where(RemoteAssetLink.remote_id == remote_id)
            return list(session.exec(query).all())

    def remote_status(self, remote: str | int | None = None) -> dict[int, dict]:
        """
        Get push status for remotes.

        Returns a dict of remote_id ->
        {"remote": Remote, "pushed": [links], "unpushed": [links]}.
        Remotes with no asset links are included, with empty lists.

        Args:
            remote: Optional filter, by remote ID or by remote name (description).
        """
        with self._database_session as session:
            query = select(Remote)
            if remote is not None:
                try:
                    query = query.where(Remote.id == int(remote))
                except ValueError:
                    query = query.where(Remote.description == remote)
            remotes = session.exec(query).unique().all()
            if not remotes:
                return {}

            by_remote = {r.id: {"remote": r, "pushed": [], "unpushed": []} for r in remotes}
            links = session.exec(
                select(RemoteAssetLink).where(RemoteAssetLink.remote_id.in_(by_remote))
            ).all()
        for link in links:
            by_remote[link.remote_id]["pushed" if link.pushed else "unpushed"].append(link)
        return by_remote

    def remote_table(self) -> Table:
        """
        Return the remote table with linked asset counts.

        Returns:
            Table: The remote table.
        """
        with self._database_session as session:
            remotes = session.exec(select(Remote)).unique().all()

            # Get link counts per remote
            links = session.exec(select(RemoteAssetLink)).all()
            pushed_counts = defaultdict(int)
            unpushed_counts = defaultdict(int)
            for link in links:
                if link.pushed:
                    pushed_counts[link.remote_id] += 1
                else:
                    unpushed_counts[link.remote_id] += 1

            return build_table(
                "Remotes",
                ["ID", "Type", "Prefix", "Description", "Push Command", "Pushed", "Unpushed"],
                [
                    (
                        str(remote.id),
                        remote.type.value,
                        remote.prefix,
                        remote.description or "No description",
                        remote.push_command or "-",
                        str(pushed_counts.get(remote.id, 0)),
                        str(unpushed_counts.get(remote.id, 0)),
                    )
                    for remote in remotes
                ],
            )

    def get_server_subscriptions(self) -> Iterable[str]:
        """
        Get the list of server subscriptions.

        Returns:
            list[str]: The list of server subscriptions.
        """
        return self.get_latest().servers

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
                select(Configuration)
                .options(selectinload(Configuration.remotes))
                .order_by(desc(Configuration.id)),
                MissingConfigDataError("No configuration row in the database"),
                unique=True,
            )
