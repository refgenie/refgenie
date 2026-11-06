"""Where refgenie reads its environment variables and builds the process-wide ``config``.

Runs once, at import: reads the ``REFGENIE_*`` variables into module constants,
creates the home directory if it is missing, constructs ``config`` from those
values, and warns if the legacy ``$REFGENIE`` variable is still set. Reload
this module (not the ``refgenie.config`` package) to re-read the environment.
"""

import os
import sys
from pathlib import Path

from refgenie.config.settings import RefgenieConfig

REFGENIE_HOME_PATH: Path = Path(os.environ.get("REFGENIE_HOME_PATH", Path.home() / ".refgenie"))
if not REFGENIE_HOME_PATH.exists():
    REFGENIE_HOME_PATH.mkdir(parents=True, exist_ok=True)

REFGENIE_LOG_LEVEL = os.environ.get("REFGENIE_LOG_LEVEL", "INFO")
REFGENIE_GENOME_FOLDER = Path(
    os.environ.get("REFGENIE_GENOME_FOLDER", REFGENIE_HOME_PATH / "genomes")
)
REFGENIE_GENOME_STAGE_FOLDER = Path(
    os.environ.get("REFGENIE_GENOME_STAGE_FOLDER", REFGENIE_HOME_PATH / "archives")
)
REFGENIE_DB_CONFIG_PATH = Path(
    os.environ.get("REFGENIE_DB_CONFIG_PATH", REFGENIE_HOME_PATH / "refgenie_db_config.yaml")
)

# Plugin switches only. Plugin *settings* live in the database
# (`refgenie plugins set`), never in environment variables.
REFGENIE_DISABLE_PLUGINS = os.environ.get("REFGENIE_DISABLE_PLUGINS", "")
REFGENIE_SERVER_PLUGINS = os.environ.get("REFGENIE_SERVER_PLUGINS", "false").lower() in (
    "1",
    "true",
    "yes",
)
# Whether a newly created configuration subscribes to the default server.
REFGENIE_SUBSCRIBE_DEFAULT_SERVER = os.environ.get(
    "REFGENIE_SUBSCRIBE_DEFAULT_SERVER", "true"
).lower() not in ("0", "false", "no")
# What runs a recipe's container image: docker, apptainer, or auto.
REFGENIE_CONTAINER_ENGINE = os.environ.get("REFGENIE_CONTAINER_ENGINE", "auto")
REFGENIE_APPTAINER_CACHEDIR = os.environ.get("REFGENIE_APPTAINER_CACHEDIR") or None

config = RefgenieConfig(
    log_level=REFGENIE_LOG_LEVEL,
    genome_folder=REFGENIE_GENOME_FOLDER,
    genome_stage_folder=REFGENIE_GENOME_STAGE_FOLDER,
    database_config_path=REFGENIE_DB_CONFIG_PATH,
    disable_plugins=REFGENIE_DISABLE_PLUGINS,
    server_plugins=REFGENIE_SERVER_PLUGINS,
    subscribe_default_server=REFGENIE_SUBSCRIBE_DEFAULT_SERVER,
    container_engine=REFGENIE_CONTAINER_ENGINE,
    apptainer_cachedir=REFGENIE_APPTAINER_CACHEDIR,
)

# Legacy refgenconf read $REFGENIE; warn if it is set without the new variable.
if "REFGENIE" in os.environ and "REFGENIE_DB_CONFIG_PATH" not in os.environ:
    _legacy = os.environ["REFGENIE"]
    print(
        "[refgenie] WARNING: $REFGENIE is set but is no longer read by refgenie 1.0.\n"
        f"           legacy value: {_legacy}\n"
        "           refgenie 1.0 uses $REFGENIE_DB_CONFIG_PATH instead.\n"
        "           See the README 'Environment variables' section for the full migration map.",
        file=sys.stderr,
    )
