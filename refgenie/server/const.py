"""Web-layer constants for the refgenie app factory (both modes)."""

import os
from importlib.metadata import version
from typing import Literal

from refgenie.const import DEFAULT_SERVER_URL

ALL_VERSIONS = {"version": version("refgenie")}

# The base URL clients use to reach this server. It comes from an environment
# variable, not from the request's Host, because behind Cloudflare and ECS
# uvicorn sees plain http.
REFGENIE_PUBLIC_URL = os.environ.get("REFGENIE_PUBLIC_URL", DEFAULT_SERVER_URL).rstrip("/")

# App modes: `server` is the public REST API; `local` is the single-user app
# behind `refgenie dash`. Unrelated to refgenie/core/mode.py's mode objects.
APP_MODE_SERVER = "server"
APP_MODE_LOCAL = "local"
AppMode = Literal["server", "local"]

#: Path prefixes owned by the API. An unmatched path under one gets a JSON 404, not index.html.
API_PATH_PREFIXES: tuple[str, ...] = (
    "/v4/",
    "/v1/",
    "/seqcol",
    "/ga4gh/",
    "/data_channel",
    "/mcp",
    "/openapi.json",
    "/docs",
    "/redoc",
    "/service-info",
    "/ping",
    "/_app/",
)

# Service-info identity shared by /service-info and the seqcol sub-service.
SERVICE_ORGANIZATION = {"name": "Refgenie", "url": "https://refgenie.org"}
SERVICE_CONTACT_URL = "https://github.com/refgenie/refgenie/issues"
SERVICE_DOCUMENTATION_URL = "https://docs.refgenie.org/refgenie/"
SERVICE_GITHUB_URL = "https://github.com/refgenie/refgenie"

# Seqcol is mounted off the root so refget's routes do not collide with refgenie's.
SEQCOL_MOUNT_PATH = "/seqcol"
SEQCOL_SERVICE_ID = "org.refgenie.seqcol"
SEQCOL_SERVICE_NAME = "Refgenie Sequence Collections"
DOWNLOAD_COUNT_DUMP_JOB_INTERVAL_SECONDS = int(
    os.environ.get("DOWNLOAD_COUNT_DUMP_JOB_INTERVAL_SECONDS", 60)
)

# Publish-catalog artifact to import at startup and daily. Unset: serve the
# SQL catalog as it is.
REFGENIE_CATALOG_URL = os.environ.get("REFGENIE_CATALOG_URL", None)
CATALOG_IMPORT_INTERVAL_SECONDS = int(os.environ.get("CATALOG_IMPORT_INTERVAL_SECONDS", 24 * 3600))
