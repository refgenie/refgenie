"""Database lifecycle: init, migrations, purge, and catalog-wide counts."""

import shutil
from pathlib import Path

from rich import print as rprint

from refgenie.config import config
from refgenie.const import DEFAULT_SERVER_URL, TARGET_ALEMBIC_VERSION
from refgenie.db import schema
from refgenie.db.engine import Engine
from refgenie.db.migrations.utils import run_sql_migrations
from refgenie.db.queries import configuration as configuration_queries
from refgenie.db.queries import summary as summary_queries
from refgenie.db.tables import Configuration
from refgenie.exceptions import DatabaseMigrationError
from refgenie.logger import logger
from refgenie.managers.base import ResourceManager
from refgenie.managers.configuration import ConfigurationManager
from refgenie.managers.sources.genomes import normalize_server_url
from refgenie.utils.prompt import Confirmer, resolve_confirmer


def _fmt_size(size_bytes: int) -> str:
    """Format bytes into human-readable size."""
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size_bytes < 1024:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024
    return f"{size_bytes:.1f} PB"


class DatabaseManager(ResourceManager):
    """Schema creation, backend init, migrations, and purge for one database."""

    def __init__(self, database_engine: Engine, configuration: ConfigurationManager):
        super().__init__(database_engine)
        self._configuration = configuration

    def summary(self) -> dict[str, int]:
        """How many genomes, asset groups and assets the catalog holds."""
        with self._database_session as session:
            return summary_queries.count_all(session)

    def species_summary(self) -> dict[str, dict[str, int]]:
        """
        Genomes, asset classes and assets per species, for species with assets.

        Returns:
            ``{species: {"genomes", "asset_classes", "assets"}}``. Genomes with
            no species name are grouped under ``"Unknown"``.
        """
        with self._database_session as session:
            rows = summary_queries.by_species(session)
        return {
            species or "Unknown": {"genomes": genomes, "asset_classes": classes, "assets": assets}
            for species, genomes, classes, assets in rows
        }

    def create_tables(self) -> None:
        """Create the database and tables; see :func:`refgenie.db.schema.create_tables`."""
        schema.create_tables(self.database_engine)

    def init(
        self,
        genome_folder: Path | None = None,
        genome_stage_folder: Path | None = None,
        config_version: int | None = None,
        subscribe_default_server: bool | None = None,
    ) -> bool:
        """
        Initialize refgenie backend and configuration.

        A new configuration is subscribed to ``DEFAULT_SERVER_URL`` so the
        first pull works without a separate ``refgenie subscribe``. Only the
        URL is recorded; nothing is fetched. An already-initialized backend
        keeps its subscriptions untouched.

        Args:
            genome_folder: The path to the genome folder, where the assets files are stored.
            genome_stage_folder: The path to the stage folder, where the genome staged assets are stored.
            config_version: Configuration row version (defaults to 1).
            subscribe_default_server: Subscribe a new configuration to
                ``DEFAULT_SERVER_URL``. None defers to
                ``config.subscribe_default_server``
                (``REFGENIE_SUBSCRIBE_DEFAULT_SERVER``, default true).

        Returns:
            bool: True if the backend is initialized, False if the Configuration
                row could not be written.
        """
        resolved_genome_folder = Path(genome_folder or config.genome_folder)
        resolved_stage_folder = (
            Path(genome_stage_folder)
            if genome_stage_folder is not None
            else (Path(config.genome_stage_folder) if config.genome_stage_folder else None)
        )

        resolved_genome_folder.mkdir(parents=True, exist_ok=True)
        logger.info(f"Genome folder ready: {resolved_genome_folder}")

        if resolved_stage_folder is not None:
            resolved_stage_folder.mkdir(parents=True, exist_ok=True)
            logger.info(f"Genome stage folder ready: {resolved_stage_folder}")

        self.create_tables()

        with self._database_session as session:
            # Idempotent: an existing Configuration row means already initialized.
            if configuration_queries.exists_any(session):
                logger.debug("refgenie backend already initialized; skipping Configuration insert")
                return True
            if subscribe_default_server is None:
                subscribe_default_server = config.subscribe_default_server
            servers = [normalize_server_url(DEFAULT_SERVER_URL)] if subscribe_default_server else []
            session.add(
                Configuration(
                    genome_folder=resolved_genome_folder.as_posix(),
                    version=config_version or 1,
                    servers=servers,
                    genome_stage_folder=(
                        resolved_stage_folder.as_posix() if resolved_stage_folder else None
                    ),
                )
            )
            try:
                session.commit()
            except Exception as e:
                logger.error(
                    f"Failed to initialize refgenie backend: {e.__class__.__name__}. "
                    "If you'd like to override current configuration and data, use 'refgenie purge'"
                )
                logger.debug(e)
                session.rollback()
                return False
        if servers:
            logger.info(f"Subscribed to the default server: {servers[0]}")
        return True

    def _current_alembic_revision(self) -> str | None:
        """The alembic revision the database is currently at, or None."""
        return schema.current_revision(self.database_engine)

    def needs_migration(self, log: bool = True) -> bool:
        """
        Check whether the database schema is behind TARGET_ALEMBIC_VERSION.

        Side effect: a database with no alembic revision is initialized first
        (``self.init()`` with the configured folders).

        Returns:
            bool: True if migrations are needed.
        """
        if (current_rev := self._current_alembic_revision()) is None:
            logger.warning(
                "Could not determine the current revision of the database schema. "
                "Running refgenie init to initialize the database"
            )
            self.init()
            current_rev = self._current_alembic_revision()
        if (target := TARGET_ALEMBIC_VERSION) is None or target == current_rev:
            logger.debug(f"Database schema is up-to-date: {current_rev}")
            return False

        if log:
            logger.warning("Database schema is outdated. Please apply migrations.")
            logger.info(
                f"Current database schema version: {current_rev}. "
                f"Required database schema version: {target}"
            )
        return True

    def migrate(self) -> None:
        """
        Migrate the database to the required version, if it is behind.

        Raises:
            DatabaseMigrationError: if the migration fails. Running on a schema
                the code does not match only fails later, at some query.
        """
        if not self.needs_migration():
            logger.debug("Database schema is up-to-date")
            return
        current = self._current_alembic_revision()
        try:
            run_sql_migrations(db_url=self.database_engine.url, revision=TARGET_ALEMBIC_VERSION)
        except Exception as e:
            raise DatabaseMigrationError(current, TARGET_ALEMBIC_VERSION, e) from e
        logger.info("Database schema updated successfully")

    def purge(self, force: bool = False, confirm: Confirmer | None = None) -> None:
        """
        Purge the refgenie backend

        Args:
            force: Whether to force the purge.
            confirm: Confirmation callback. Defaults to a refusal unless the CLI
                has enabled interactive prompts; see `refgenie.utils.prompt`.
        """
        if not schema.has_schema(self.database_engine):
            logger.info(
                "No database tables found. Nothing to purge. "
                "You may want to check if the database is set up correctly. "
                "If you're sure it is, purge the refgenie-managed files manually."
            )
            return

        if not force:
            # Read with raw SQL, so a schema mismatch still gets a preview.
            genome_folder, genomes = schema.purge_preview(self.database_engine)
            genome_folder = genome_folder or config.genome_folder

            rprint(f"\n  Database: {self.database_engine.url}")
            rprint(f"  Genome folder: {genome_folder}")
            if genomes:
                rprint("  Genomes:")
                for name, asset_count, total_size in genomes:
                    size = _fmt_size(total_size)
                    rprint(f"    {name}: {asset_count} assets ({size})")
            else:
                rprint("  Genomes: (none)")
            rprint()

            if not resolve_confirmer(confirm, default=False)("Purge all of this?"):
                logger.info("Aborted by a user")
                return

        try:
            genome_folder = Path(self._configuration.get_latest().genome_folder)
        except Exception:
            genome_folder = config.genome_folder

        # Schema first, files second: an interrupted purge leaves inert files,
        # never a catalog pointing at deleted data.
        schema.drop_tables(self.database_engine)
        logger.info(f"Purged refgenie backend: '{self.database_engine.url}'")
        if genome_folder.exists():
            shutil.rmtree(genome_folder, ignore_errors=True)
            logger.info(f"Purged genome folder: '{genome_folder}'")
