"""Database lifecycle: config, engine creation, init, migrations, and purge."""

import shutil
from pathlib import Path

from alembic.migration import MigrationContext
from rich import print as rprint
from sqlalchemy.engine import Engine as SqlalchemyDatabaseEngine
from sqlmodel import SQLModel, create_engine, select

from refgenie.config import config
from refgenie.config.db import (
    DatabaseConfig,
    create_default_db_config,
    dump_db_config,
    load_db_config,
)
from refgenie.const import TARGET_ALEMBIC_VERSION
from refgenie.db.migrations.utils import run_sql_migrations
from refgenie.db.tables import AlembicVersion, Configuration
from refgenie.logger import logger
from refgenie.utils.prompt import Confirmer, resolve_confirmer


def _fmt_size(size_bytes: int) -> str:
    """Format bytes into human-readable size."""
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size_bytes < 1024:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024
    return f"{size_bytes:.1f} PB"


class DatabaseLifecycleMixin:
    """
    Database config/engine construction, backend init, migrations, and purge.

    Mixed into :class:`refgenie.core.facade.Refgenie`; relies on the facade's
    ``database_engine``, ``_database_session``, and manager properties.
    """

    @staticmethod
    def get_database_config(config_path: str | Path | None = None) -> DatabaseConfig:
        """
        Get the database configuration.

        Args:
            config_path: The path to the database configuration file.

        Returns:
            DatabaseConfig: The database configuration.
        """
        if isinstance(config_path, str):
            config_path = Path(config_path)
        cp = config_path or config.database_config_path
        if isinstance(cp, str):
            cp = Path(cp)
        if not cp.exists():
            logger.debug(
                f"Database configuration file not found at {cp}. "
                "Creating a new one with default values."
            )
            cp.parent.mkdir(parents=True, exist_ok=True)
            cp.write_text(dump_db_config(create_default_db_config()))
            logger.info(
                f"Database configuration file created at {cp}. "
                "Please edit it to set the database connection parameters."
            )
        return load_db_config(cp)

    @classmethod
    def get_default_database_engine(
        cls, database_config_path: str | Path | None = None, echo: bool = False
    ) -> SqlalchemyDatabaseEngine:
        """
        Get the default database engine.

        Args:
            database_config_path: The path to the database configuration file, if not provided,
                the default configuration is used.
            echo: Whether to echo SQL statements.

        Returns:
            Engine: The database engine.
        """
        database_config = cls.get_database_config(config_path=database_config_path)
        # pre_ping + recycle guard against silently-dead pooled connections.
        # Observed on the deployed server (Postgres on RDS): after an idle
        # stretch, the first requests hung ~2 minutes in TCP retransmission on
        # reaped connections, then all completed at once. A liveness check per
        # checkout is trivial next to that; both options are no-ops for SQLite.
        return create_engine(
            url=database_config.url, echo=echo, pool_pre_ping=True, pool_recycle=300
        )

    def _create_db_and_tables(self):
        """
        Create the database and tables.

        New databases get the full current schema from ``create_all`` and are
        stamped at the migration head. Existing databases are brought forward
        by alembic (``check_for_db_migrations`` / ``migrate_db``), which is the
        only schema-evolution path.
        """
        SQLModel.metadata.create_all(self.database_engine)
        logger.info(f"Initialized refgenie backend: '{self._database_engine.url}'")
        with self._database_session as session:
            if session.exec(select(AlembicVersion)).first() is not None:
                return
            session.add(AlembicVersion(version_num=TARGET_ALEMBIC_VERSION))
            session.commit()
            logger.debug(f"Set alembic version to {TARGET_ALEMBIC_VERSION}")

    def init_backend(
        self,
        genome_folder: Path,
        genome_stage_folder: Path | None = None,
        config_version: int | None = None,
    ) -> bool:
        """
        Initialize refgenie backend and configuration.

        Args:
            genome_folder: The path to the genome folder, where the assets files are stored.
            genome_stage_folder: The path to the stage folder, where the genome staged assets are stored.
            config_version: Configuration row version (defaults to 1).
        """
        self._create_db_and_tables()

        # set the genome folder
        with self._database_session as session:
            # Idempotent init: the persistent build catalog re-runs init nightly,
            # and Configuration.version is unique, so re-inserting a row would
            # raise IntegrityError. If a Configuration row already exists the
            # backend is initialized; skip the insert.
            if session.exec(select(Configuration)).first() is not None:
                logger.debug("refgenie backend already initialized; skipping Configuration insert")
                return True
            session.add(
                Configuration(
                    genome_folder=genome_folder.as_posix(),
                    version=config_version or 1,
                    genome_stage_folder=(
                        genome_stage_folder.as_posix() if genome_stage_folder else None
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
        return True

    def init(
        self,
        genome_folder: Path | None = None,
        genome_stage_folder: Path | None = None,
        config_version: int | None = None,
    ):
        """
        Initialize refgenie backend and configuration.

        Args:
            genome_folder: The path to the genome folder, where the assets files are stored.
            genome_stage_folder: The path to the stage folder, where the genome staged assets are stored.
            config_version: Configuration row version (defaults to 1).
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

        self.init_backend(
            genome_folder=resolved_genome_folder,
            genome_stage_folder=resolved_stage_folder,
            config_version=config_version,
        )

    def _current_alembic_revision(self) -> str | None:
        """
        Read the alembic revision the database is currently at.

        Runs on every Refgenie construction; keep the connection
        context-managed so it returns to the pool.

        Returns:
            str | None: The current revision, or None if the database has none.
        """
        with self.database_engine.connect() as connection:
            return MigrationContext.configure(connection).get_current_revision()

    def check_for_db_migrations(self, log: bool = True) -> bool:
        """
        Check whether the database schema is behind TARGET_ALEMBIC_VERSION,
        initializing the database first if it has no revision.

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

    def migrate_db(self):
        """
        Migrate the database to the required version.
        """
        if not self.check_for_db_migrations():
            logger.info("Database schema is up-to-date")
            return
        try:
            run_sql_migrations(db_url=self.database_engine.url, revision=TARGET_ALEMBIC_VERSION)
        except Exception as e:
            logger.exception(
                f"Failed to apply database migrations: {e}. "
                "Refgenie configuration may be corrupted. Please contact the developers."
            )
        else:
            logger.info("Database schema updated successfully")

    def check_table_exists(self, table_name: str | None = None) -> bool:
        """
        Check if the tables exist.

        Returns:
            bool: Whether the tables exist.
        """
        from sqlalchemy import inspect

        return inspect(self.database_engine).has_table(table_name or "configuration")

    def purge(self, force: bool = False, confirm: Confirmer | None = None):
        """
        Purge the refgenie backend

        Args:
            force: Whether to force the purge.
            confirm: Confirmation callback. Defaults to a refusal unless the CLI
                has enabled interactive prompts; see `refgenie.utils.prompt`.
        """
        if not self.check_table_exists():
            logger.info(
                "No database tables found. Nothing to purge. "
                "You may want to check if the database is set up correctly. "
                "If you're sure it is, purge the refgenie-managed files manually."
            )
            return

        if not force:
            from sqlalchemy import text

            # Get genome folder from DB (resilient to schema mismatches)
            genome_folder = config.genome_folder
            try:
                with self._database_engine.connect() as conn:
                    row = conn.execute(
                        text("SELECT genome_folder FROM configuration ORDER BY id DESC LIMIT 1")
                    ).first()
                    if row:
                        genome_folder = Path(row[0])
            except Exception:
                pass

            # Per-genome summary
            genomes = []
            try:
                with self._database_engine.connect() as conn:
                    genomes = conn.execute(
                        text(
                            "SELECT a.name, COUNT(DISTINCT asset.digest), COALESCE(SUM(asset.size), 0) "
                            "FROM genome g "
                            "LEFT JOIN alias a ON a.genome_digest = g.digest "
                            "LEFT JOIN assetgroup ag ON ag.genome_digest = g.digest "
                            "LEFT JOIN asset ON asset.asset_group_id = ag.id "
                            "GROUP BY g.digest, a.name"
                        )
                    ).fetchall()
            except Exception:
                pass

            rprint(f"\n  Database: {self._database_engine.url}")
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

        # Get genome folder before dropping tables
        try:
            genome_folder = self.genome_folder
        except Exception:
            genome_folder = config.genome_folder

        # Schema first, files second. An interrupted purge should leave data with
        # no catalog -- inert bytes a later purge or a manual rm clears -- rather
        # than a catalog describing data that is already gone, which every read
        # path would believe.
        SQLModel.metadata.drop_all(self._database_engine)
        logger.info(f"Purged refgenie backend: '{self._database_engine.url}'")
        if genome_folder.exists():
            shutil.rmtree(genome_folder, ignore_errors=True)
            logger.info(f"Purged genome folder: '{genome_folder}'")
