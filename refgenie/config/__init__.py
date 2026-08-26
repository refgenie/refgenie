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

# --- Localhost bridge (refgenie dash <-> the public SPA) ---------------------
#
# Whether and how much a browser page on an allowlisted public origin may talk
# to a locally running `refgenie dash`. Consumed by
# refgenie/server/local/security.py (LocalSecuritySettings reads the same env
# vars at app construction; these constants are its defaults and the CLI's
# view of the configuration).
#
#   off  -- no cross-origin access at all (same-origin dash UI still works)
#   read -- allowlisted origins may read (/ping, /v4, /v1/jobs); no actions
#   full -- additionally allows POST /v1/actions/pull cross-origin
REFGENIE_BRIDGE_MODE = os.environ.get("REFGENIE_BRIDGE_MODE", "read")
#: Comma-separated exact origins. The deployed SPA lives at the apex
#: https://refgenie.org; https://ui.refgenie.org is its former home and still
#: 301-redirects there. The docs site is https://docs.refgenie.org and is
#: deliberately NOT allowlisted.
REFGENIE_BRIDGE_ORIGINS = os.environ.get(
    "REFGENIE_BRIDGE_ORIGINS", "https://refgenie.org,https://ui.refgenie.org"
)
#: Optional origin regex. Exists for local SPA development and Cloudflare
#: preview origins ONLY. A careless regex here (e.g. ".*") is the single
#: easiest way to undo every protection the bridge design provides.
REFGENIE_BRIDGE_ORIGIN_REGEX = os.environ.get("REFGENIE_BRIDGE_ORIGIN_REGEX", "")
#: Whether /ping may reveal filesystem paths. Default false: a home-directory
#: path leaks the OS username to every allowlisted origin.
REFGENIE_BRIDGE_EXPOSE_PATHS = os.environ.get("REFGENIE_BRIDGE_EXPOSE_PATHS", "false").lower() in (
    "1",
    "true",
    "yes",
)

config = RefgenieConfig(
    log_level=REFGENIE_LOG_LEVEL,
    genome_folder=REFGENIE_GENOME_FOLDER,
    genome_stage_folder=REFGENIE_GENOME_STAGE_FOLDER,
    database_config_path=REFGENIE_DB_CONFIG_PATH,
    bridge_mode=REFGENIE_BRIDGE_MODE,
    bridge_origins=REFGENIE_BRIDGE_ORIGINS,
    bridge_origin_regex=REFGENIE_BRIDGE_ORIGIN_REGEX,
    bridge_expose_paths=REFGENIE_BRIDGE_EXPOSE_PATHS,
)

# Legacy env-var migration warning.
# Legacy refgenconf read $REFGENIE; refgenie1 uses $REFGENIE_DB_CONFIG_PATH (and friends).
# If a user still has $REFGENIE set but has not adopted the new variable, warn
# loudly rather than silently ignoring it.
if "REFGENIE" in os.environ and "REFGENIE_DB_CONFIG_PATH" not in os.environ:
    _legacy = os.environ["REFGENIE"]
    print(
        "[refgenie] WARNING: $REFGENIE is set but is no longer read by refgenie1.\n"
        f"           legacy value: {_legacy}\n"
        "           refgenie1 uses $REFGENIE_DB_CONFIG_PATH instead.\n"
        "           See the README 'Environment variables' section for the full migration map.",
        file=sys.stderr,
    )
