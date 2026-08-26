from string import ascii_letters, digits

DEFAULT_PULL_SIZE_CUTOFF_GB = 10
DEFAULT_PAGE_SIZE = 100
MAX_PAGE_SIZE = 1000

ASSET_NAME_CHAR_WHITELIST = set(ascii_letters + digits + "-._~")
ASSET_REGISTRY_PATH_COMPONENT_CHAR_BLACKLIST = {" ", ":", "/"}

# Top-level tree under genome_folder holding per-build bookkeeping (logs,
# commands, profile, stats, completion flags). Deliberately NOT inside the
# asset directory: build records describe an event, not asset content.
BUILDS_DIR = "builds"

# Top-level tree under genome_folder holding the name-addressed view of the
# digest-addressed data/ tree: alias/<alias>/<group>/<asset>/ of symlinks.
ALIAS_DIR = "alias"

API_VERSION = "v4"

#: Server offered to a fresh instance with no subscriptions (facade.pull's
#: subscribe prompt). The v4 server refgenie1 actually speaks to -- not the
#: legacy v3 refgenomes.databio.org, which 404s on /service-info.
DEFAULT_SERVER_URL = "https://api.refgenie.org"

# CWE-321: Hardcoded key retained for demo/trial convenience.
# Local-only DB; attacker with filesystem access already has everything.
# Override with REFGENIE_ENCRYPTION_KEY env var for production.
DEFAULT_ENCRYPTION_KEY = "MkNSaEhqMU9SS2dzSFBIT3g3bkdwcWZ3Y2l0akhqQjFZd2t6RjladEJWZz0="

CURRENT_CONFIG_VERSION = 1

#: Head revision in refgenie/db/migrations/versions. New databases are stamped
#: with this; existing ones are migrated up to it. Update whenever a migration
#: is added -- it must always name the head, or later migrations never run.
TARGET_ALEMBIC_VERSION = "b2f1c3d4e5a6"


def __getattr__(name: str):
    """Compute version-derived constants lazily.

    ``importlib.metadata.version()`` scans the environment's installed
    distributions, which is the only non-trivial work this module would
    otherwise do at import time. Deferring it keeps ``import refgenie.const``
    a sub-millisecond import of pure literals.
    """
    if name in ("REFGENIE_VERSION", "USER_AGENT", "HTTP_HEADERS"):
        from importlib.metadata import version

        refgenie_version = version("refgenie")
        # A non-default User-Agent is required to get through Cloudflare (which
        # fronts api.refgenie.org and blocks the default Python-urllib /
        # python-httpx UAs with a 403). Naming it after the package is honest
        # and debuggable.
        user_agent = f"refgenie/{refgenie_version}"
        values = {
            "REFGENIE_VERSION": refgenie_version,
            "USER_AGENT": user_agent,
            "HTTP_HEADERS": {"User-Agent": user_agent},
        }
        globals().update(values)
        return values[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
