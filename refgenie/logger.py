import logging

from rich.logging import RichHandler

from refgenie.config import config


def get_logger(name: str = "refgenie", level: str = config.log_level.value) -> logging.Logger:
    """
    Get a logger with the given name.
    """
    logging.getLogger().handlers.clear()
    logging.basicConfig(
        level=level,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(show_time=False, rich_tracebacks=True)],
    )
    logging.getLogger("alembic").setLevel(logging.WARNING)
    return logging.getLogger(name)


logger = get_logger()
