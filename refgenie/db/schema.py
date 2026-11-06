"""Creating, inspecting and dropping the refgenie schema."""

from pathlib import Path

from alembic.migration import MigrationContext
from sqlalchemy import inspect, text
from sqlalchemy.engine import Engine
from sqlmodel import SQLModel, select

from refgenie.const import TARGET_ALEMBIC_VERSION
from refgenie.db.session import database_session
from refgenie.db.tables import AlembicVersion
from refgenie.logger import logger


def create_tables(engine: Engine) -> None:
    """
    Create any missing tables, and stamp a new database at the migration head.

    Existing databases are upgraded only by alembic.
    """
    SQLModel.metadata.create_all(engine)
    logger.info(f"Initialized refgenie backend: '{engine.url}'")
    with database_session(engine) as session:
        if session.exec(select(AlembicVersion)).first() is not None:
            return
        session.add(AlembicVersion(version_num=TARGET_ALEMBIC_VERSION))
        session.commit()
        logger.debug(f"Set alembic version to {TARGET_ALEMBIC_VERSION}")


def drop_tables(engine: Engine) -> None:
    """Drop every refgenie table."""
    SQLModel.metadata.drop_all(engine)


def has_schema(engine: Engine) -> bool:
    """Whether the refgenie schema (the ``configuration`` table) exists."""
    return inspect(engine).has_table("configuration")


def current_revision(engine: Engine) -> str | None:
    """The alembic revision the database is at, or None if it has none."""
    with engine.connect() as connection:
        return MigrationContext.configure(connection).get_current_revision()


def purge_preview(engine: Engine) -> tuple[Path | None, list[tuple[str, int, int]]]:
    """
    What a purge would remove: the genome folder and a per-genome summary.

    Raw SQL, so a database whose schema does not match the code can still be
    described. A half that fails comes back empty rather than raising.

    Returns:
        The genome folder recorded in the newest configuration row, and
        ``(alias, asset count, total size)`` rows.
    """
    genome_folder = None
    try:
        with engine.connect() as conn:
            row = conn.execute(
                text("SELECT genome_folder FROM configuration ORDER BY id DESC LIMIT 1")
            ).first()
            if row:
                genome_folder = Path(row[0])
    except Exception:
        pass

    genomes = []
    try:
        with engine.connect() as conn:
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
    return genome_folder, [tuple(row) for row in genomes]
