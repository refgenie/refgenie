from pathlib import Path

from alembic import command
from alembic.config import Config
from sqlalchemy.engine.url import URL

from refgenie.logger import logger


def run_sql_migrations(db_url: URL, revision: str = "head") -> None:
    # Config() is deliberately empty: alembic.ini is a development-tree file and
    # is not shipped in the wheel, so every option is set explicitly below.
    config = Config()
    # The migrations directory is where this utils.py file is located
    migrations_dir = Path(__file__).parent
    config.set_main_option("script_location", str(migrations_dir))
    # set_main_option stores into a ConfigParser, whose BasicInterpolation reads
    # `%` as interpolation syntax. A URL whose password percent-encodes to %XX
    # (special chars) then fails with "invalid interpolation syntax". Escape % as
    # %% so ConfigParser restores the literal URL when env.py reads it back.
    url_str = db_url.render_as_string(hide_password=False).replace("%", "%%")
    config.set_main_option("sqlalchemy.url", url_str)

    # Set other necessary options from alembic.ini manually
    config.set_main_option("prepend_sys_path", ".")
    config.set_main_option("version_path_separator", "os")

    logger.debug(f"Running migrations from {migrations_dir} to revision '{revision}'")

    command.upgrade(config=config, revision=revision)
    logger.info("Migrations complete.")
