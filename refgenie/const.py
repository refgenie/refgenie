from string import ascii_letters, digits

DEFAULT_PULL_SIZE_CUTOFF_GB = 10
DEFAULT_PAGE_SIZE = 100
MAX_PAGE_SIZE = 1000
# PostgreSQL INTEGER and BIGINT (OFFSET) limits; query params are capped at
# these so an oversized number is a 422, not a 500.
SQL_INT_MAX = 2**31 - 1
SQL_BIGINT_MAX = 2**63 - 1

ASSET_NAME_CHAR_WHITELIST = set(ascii_letters + digits + "-._~")
ASSET_REGISTRY_PATH_COMPONENT_CHAR_BLACKLIST = {" ", ":", "/"}

# Tree under genome_folder for per-build records (logs, commands, stats, flags),
# kept outside the asset directory because they are not asset content.
BUILDS_DIR = "builds"

# Top-level tree under genome_folder holding the name-addressed view of the
# digest-addressed data/ tree: alias/<alias>/<group>/<asset>/ of symlinks.
ALIAS_DIR = "alias"

API_VERSION = "v4"

#: Server offered to a fresh instance with no subscriptions (the v4 API).
DEFAULT_SERVER_URL = "https://api.refgenie.org"

# CWE-321: hardcoded key for a local-only DB, where filesystem access already
# exposes everything. Override with REFGENIE_ENCRYPTION_KEY.
DEFAULT_ENCRYPTION_KEY = "MkNSaEhqMU9SS2dzSFBIT3g3bkdwcWZ3Y2l0akhqQjFZd2t6RjladEJWZz0="

CURRENT_CONFIG_VERSION = 1

#: Head revision in refgenie/db/migrations/versions. New databases are stamped
#: with it and old ones migrated to it, so update it with every new migration.
TARGET_ALEMBIC_VERSION = "05b8e0386023"


def __getattr__(name: str):
    """Compute version-derived constants lazily, since reading metadata is slow."""
    if name in ("REFGENIE_VERSION", "USER_AGENT", "HTTP_HEADERS"):
        from importlib.metadata import version

        refgenie_version = version("refgenie")
        # Cloudflare in front of api.refgenie.org blocks default client UAs.
        user_agent = f"refgenie/{refgenie_version}"
        values = {
            "REFGENIE_VERSION": refgenie_version,
            "USER_AGENT": user_agent,
            "HTTP_HEADERS": {"User-Agent": user_agent},
        }
        globals().update(values)
        return values[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
