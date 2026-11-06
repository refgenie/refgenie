from pathlib import Path

from alembic import command
from alembic.config import Config
from sqlalchemy.engine.url import URL

from refgenie.logger import logger


def run_sql_migrations(db_url: URL, revision: str = "head") -> None:
    # alembic.ini is not shipped in the wheel, so every option is set here.
    config = Config()
    migrations_dir = Path(__file__).parent
    config.set_main_option("script_location", str(migrations_dir))
    # ConfigParser treats `%` as interpolation, so escape a percent-encoded password.
    url_str = db_url.render_as_string(hide_password=False).replace("%", "%%")
    config.set_main_option("sqlalchemy.url", url_str)

    config.set_main_option("prepend_sys_path", ".")
    config.set_main_option("path_separator", "os")

    logger.debug(f"Running migrations from {migrations_dir} to revision '{revision}'")

    command.upgrade(config=config, revision=revision)
    logger.info("Migrations complete.")
