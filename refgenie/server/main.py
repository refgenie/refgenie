try:
    # the following imports are only included in extras 'dash'
    import uvicorn
    from fastapi import Depends, FastAPI
    from fastapi.middleware.cors import CORSMiddleware
    from fastapi.responses import PlainTextResponse
except ImportError as e:
    raise ImportError(
        "The 'dash' extras are not installed. Please install refgenie with the 'dash' extras to use the web interface."
    ) from e
import logging
from contextlib import asynccontextmanager
from pathlib import Path

from refgenie.server.dependencies import get_refgenie as _get_refgenie
from refgenie.server.dependencies import reject_nul_bytes
from refgenie.server.routers.catalog import router as refgenie_core_router
from refgenie.core import Refgenie

from refgenie.const import API_VERSION
from refgenie.server.const import (
    ALL_VERSIONS,
    APP_MODE_LOCAL,
    APP_MODE_SERVER,
    AppMode,
    CATALOG_IMPORT_INTERVAL_SECONDS,
    DOWNLOAD_COUNT_DUMP_JOB_INTERVAL_SECONDS,
    REFGENIE_CATALOG_URL,
    SEQCOL_MOUNT_PATH,
    SEQCOL_SERVICE_ID,
    SEQCOL_SERVICE_NAME,
    SERVICE_CONTACT_URL,
    SERVICE_DOCUMENTATION_URL,
    SERVICE_GITHUB_URL,
    SERVICE_ORGANIZATION,
)
from refgenie.server.seqcol import SeqcolService, SeqcolServing
from refgenie.server.spa import mount_spa, read_build_info, resolve_web_dist
from refgenie.server.stats import (
    EndpointCollector,
    EndpointHitCounterMiddleware,
    handle_endpoint_collector_hits,
)

logger = logging.getLogger(__name__)


SERVER_TAGS_METADATA = [
    {
        "name": "Version 4",
        "description": "Version 4 of the Refgenieserver API, which includes support for pangenomes and sequences.",
    },
    {
        "name": "Default",
        "description": "Default endpoints for the Refgenieserver API, which includes basic genome asset management.",
    },
    {
        "name": "Data Channel",
        "description": "Endpoints for the Data Channel, which provides access to genome assets in a data channel format.",
    },
    {
        "name": "GA4GH DRS",
        "description": "Endpoints for the GA4GH DRS (Data Repository Service), which provides access to genomic data in a standardized format. More details: [GA4GH DRS standard](https://ga4gh.github.io/data-repository-service-schemas/)",
    },
]

LOCAL_TAGS_METADATA = [
    {
        "name": "Version 4",
        "description": "The refgenie JSON API: genomes, assets, recipes and asset classes.",
    },
    {
        "name": "Actions",
        "description": "State-changing operations on this local refgenie instance (pull, build, delete, aliases, subscriptions).",
    },
    {
        "name": "Remote",
        "description": "Browsing the catalogs of the servers this instance subscribes to.",
    },
    {
        "name": "Local files",
        "description": "Viewing and downloading an asset's files from this machine's disk (same-origin only).",
    },
]


def _capabilities(mode: AppMode) -> dict[str, bool]:
    """The capability flags for ``mode``.

    Server mode only serves data; local mode manages one user's assets.
    ``recipes_write`` and ``asset_classes_write`` are always false: they take a
    server-side path or URL, so they stay CLI-only.
    """
    is_local = mode == APP_MODE_LOCAL
    return {
        "pull": is_local,
        "build": is_local,
        "delete": is_local,
        "aliases_write": is_local,
        "subscriptions": is_local,
        "data_channels": is_local,
        "recipes_write": False,
        "asset_classes_write": False,
        "remote_browse": is_local,
        "genome_init": is_local,
        "jobs": is_local,
        "jobs_cancel": is_local,
        "downloads": not is_local,
        "archives": not is_local,
        "seqcol": not is_local,
        "drs": not is_local,
        "local_files": is_local,
    }


def _configure_cors(app: FastAPI, mode: AppMode) -> None:
    """Install open, credential-less CORS for the read-only ``server`` mode.

    Local-mode CORS is set by ``refgenie.server.local.security.install_local_security``.
    """
    if mode != APP_MODE_SERVER:
        return
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],
        allow_credentials=False,
        allow_methods=["GET", "HEAD", "OPTIONS"],
        allow_headers=["*"],
    )


def _service_info_doc(
    mode: AppMode,
    root_path: str,
    build_info: dict | None,
    capabilities: dict[str, bool],
    seqcol_services: tuple[SeqcolService, ...],
) -> dict:
    """The ``/service-info`` document.

    Refgenie-specific fields, including ``capabilities``, nest under
    ``refgenie``. In server mode, ``seqcol.services`` lists every seqcol mount
    and the flat ``seqcol`` keys describe the default one. Local mode has no
    ``seqcol`` block. Built per request, since the store set can change.
    """
    is_local = mode == APP_MODE_LOCAL
    document = {
        "id": "org.refgenie.api",
        "name": "Refgenie",
        "type": {"group": "org.refgenie", "artifact": "refgenie", "version": API_VERSION},
        "description": (
            "A web interface and RESTful API for reference genome assets, "
            "with a GA4GH sequence collections service"
        ),
        "organization": SERVICE_ORGANIZATION,
        "contactUrl": SERVICE_CONTACT_URL,
        "documentationUrl": SERVICE_DOCUMENTATION_URL,
        "environment": "production",
        "version": ALL_VERSIONS["version"],
        "refgenie": {
            "mode": mode,
            "api_base": "/v4",
            "root_path": root_path,
            "refgenie_version": ALL_VERSIONS["version"],
            "service_name": "refgenie local dashboard" if is_local else "refgenie server",
            "web_ui": build_info,
            "capabilities": capabilities,
            "links": {
                "docs": SERVICE_DOCUMENTATION_URL,
                "github": SERVICE_GITHUB_URL,
                "openapi": "/openapi.json",
            },
        },
    }
    if is_local:
        return document
    services = [
        {
            "id": m.id,
            "name": m.name,
            "store": m.store,
            "url": m.url,
            "service_info": f"{m.url}/service-info",
            "refget_store": {"enabled": True, "url": m.store_url},
        }
        for m in seqcol_services
    ] or [
        {
            "id": SEQCOL_SERVICE_ID,
            "name": SEQCOL_SERVICE_NAME,
            "url": SEQCOL_MOUNT_PATH,
            "service_info": f"{SEQCOL_MOUNT_PATH}/service-info",
            "refget_store": {"enabled": False},
        }
    ]
    document["seqcol"] = {**services[0], "services": services}
    return document


def _mount_mcp(app: FastAPI, rg: Refgenie) -> None:
    """Serve the MCP server at ``/mcp``."""
    from mcp.server.transport_security import TransportSecuritySettings
    from refgenie.mcp.tools import mcp as mcp_server, set_refgenie as set_mcp_refgenie

    set_mcp_refgenie(rg, serves_downloads=True)
    # The default DNS-rebinding guard rejects every non-loopback Host; this API is public.
    mcp_app = mcp_server.streamable_http_app(
        transport_security=TransportSecuritySettings(enable_dns_rebinding_protection=False)
    )
    # Captured before another create_app() replaces the session manager.
    app.state.mcp_session_manager = mcp_server.session_manager
    # Routes added directly: mounting at /mcp would serve /mcp/mcp.
    app.router.routes.extend(mcp_app.routes)


def create_app(
    mode: AppMode = APP_MODE_SERVER,
    root_path: str = "",
    refgenie_instance: "Refgenie | None" = None,
    web_dist: "Path | None" = None,
) -> FastAPI:
    """Create and configure a Refgenie FastAPI application.

    * ``server`` -- the public API: v4, GA4GH DRS, data channels, seqcol,
      ``/mcp``, download counting and the background scheduler.
    * ``local`` -- ``refgenie dash``: v4 over the user's own refgenie, plus the
      command surface under ``/v1``.

    Both serve the same web UI, which reads its mode from ``/service-info``.
    This ``mode`` is independent of ``rg.mode``.

    Args:
        mode: ``"server"`` (default) or ``"local"``.
        root_path: URL prefix for sub-path mounting (e.g. "/refgenie").
            Passed to FastAPI(root_path=...) and used to rewrite the SPA's
            ``<base href>`` so a sub-path deployment needs no frontend rebuild.
        refgenie_instance: Optional pre-built Refgenie instance. If provided,
            it is used instead of creating a new one (useful for testing).
        web_dist: Optional path to a built web UI bundle, overriding the
            packaged ``refgenie/server/webui/``. Tests pass a tmp dir here.

    Returns:
        A fully configured FastAPI app ready to serve or mount as a sub-app.
    """
    if mode not in (APP_MODE_SERVER, APP_MODE_LOCAL):
        raise ValueError(f"Unknown app mode: {mode!r}")
    is_local = mode == APP_MODE_LOCAL

    if refgenie_instance is not None:
        rg = refgenie_instance
    elif is_local:
        rg = _get_refgenie()
    else:
        rg = Refgenie(server_mode=True)

    def _import_catalog():
        """Refresh the SQL catalog from the published artifact; a failure is logged, not fatal."""
        from refgenie.publish_catalog import import_publish_catalog

        try:
            import_publish_catalog(rg.database_engine, REFGENIE_CATALOG_URL)
        except Exception:
            logger.exception(f"Publish-catalog import from {REFGENIE_CATALOG_URL} failed")

    # The catalog carries the store list, so import it before anything reads
    # the stores: a fresh container learns its stores here.
    if not is_local and REFGENIE_CATALOG_URL:
        rg.database.create_tables()
        _import_catalog()

    # Seqcol follows the `store` registry: one sub-app per enabled store, or
    # DB-backed seqcol with none. The lifespan loads the stores.
    seqcol = None if is_local else SeqcolServing(rg)

    def _scheduled_catalog_import():
        """Import the catalog, then serve whatever store set it now lists."""
        _import_catalog()
        try:
            seqcol.refresh()
        except Exception:
            logger.exception("Could not load the new store set; still serving the old one")

    @asynccontextmanager
    async def local_lifespan(app: FastAPI):
        """Local mode: create the user's tables and run the JobManager while serving."""
        from refgenie.server.jobs import JobManager

        rg.database.create_tables()
        app.state.job_manager = JobManager(rg)
        try:
            yield
        finally:
            app.state.job_manager.shutdown(wait=False)
            app.state.job_manager = None

    @asynccontextmanager
    async def server_lifespan(app: FastAPI):
        """Server mode: tables, stores, seqcol, scheduler, MCP.

        The first catalog import already ran in ``create_app``. Each store is
        loaded once, read-only, and shared by the API and its seqcol sub-app;
        each scheduled import reloads them if the store set changed.
        """
        # A `server`-extra dependency; `refgenie dash` installs lack it.
        from apscheduler.schedulers.background import BackgroundScheduler

        import time

        startup = time.monotonic()
        print("Creating database and tables")
        rg.database.create_tables()
        # Load every store now; a store that fails to load stops startup.
        seqcol.refresh()
        app.state.endpoint_hits_collector = EndpointCollector()
        scheduler = BackgroundScheduler()
        scheduler.add_job(
            handle_endpoint_collector_hits,
            "interval",
            seconds=DOWNLOAD_COUNT_DUMP_JOB_INTERVAL_SECONDS,
            args=[app, rg],
        )
        if REFGENIE_CATALOG_URL:
            # Idempotent upsert, so the per-worker scheduler is harmless.
            scheduler.add_job(
                _scheduled_catalog_import,
                "interval",
                seconds=CATALOG_IMPORT_INTERVAL_SECONDS,
            )
        scheduler.start()
        logger.info(f"Server startup finished in {time.monotonic() - startup:.1f}s")
        # The /mcp sub-app's lifespan never runs, so start its session manager here.
        async with app.state.mcp_session_manager.run():
            yield
        scheduler.shutdown()

    app = FastAPI(
        title="Refgenie" if is_local else "Refgenieserver REST API",
        description=(
            "the local management API for one refgenie installation"
            if is_local
            else "a web interface and RESTful API for reference genome assets"
        ),
        version=ALL_VERSIONS["version"],
        lifespan=local_lifespan if is_local else server_lifespan,
        openapi_tags=LOCAL_TAGS_METADATA if is_local else SERVER_TAGS_METADATA,
        swagger_ui_parameters={"docExpansion": "none", "defaultModelsExpandDepth": 0},
        root_path=root_path,
        dependencies=[Depends(reject_nul_bytes)],
    )
    app.state.mode = mode

    # Registration order: middleware, dependency overrides, routers, mounts, and
    # the SPA catch-all last. Anything registered after mount_spa() is dead.

    _configure_cors(app, mode)

    # Routers read this app's Refgenie instance.
    app.dependency_overrides[_get_refgenie] = lambda: rg

    if not is_local:
        app.add_middleware(EndpointHitCounterMiddleware)

    if is_local:
        from refgenie.server.errors import install_error_handlers
        from refgenie.server.local.actions import router as actions_router
        from refgenie.server.local.security import install_local_security

        install_local_security(app)
        install_error_handlers(app)

        # Local only: some request models accept server-local filesystem paths.
        app.include_router(actions_router, prefix="/v1")

        from refgenie.server.routers.remote import router as remote_router

        app.include_router(remote_router, prefix="/v1")

        from refgenie.server.jobs import jobs_router

        app.include_router(jobs_router, prefix="/v1")

        # Read-only asset files on disk, same-origin only.
        from refgenie.server.local.files import router as local_files_router

        app.include_router(local_files_router, prefix="/v1")

    app.include_router(refgenie_core_router, prefix="/v4", tags=["Version 4"])

    if not is_local:
        # Import routers after the app is defined to avoid circular imports
        from refgenie.server.routers.version4 import router as version4

        app.include_router(version4, prefix="/v4", tags=["Version 4"])

        from refgenie.server.routers.ga4gh_drs import router as ga4gh_drs

        app.include_router(ga4gh_drs, prefix="/ga4gh/drs/v1", tags=["GA4GH DRS"])

        from refgenie.server.routers.data_channel import router as data_channel

        app.include_router(data_channel, prefix="/v4/data_channel", tags=["Data Channel"])
        app.include_router(data_channel, prefix="/data_channel", tags=["Data Channel"])

    if not is_local:
        # One mount for every store; it dispatches per request.
        app.mount(SEQCOL_MOUNT_PATH, seqcol)
        app.state.seqcol = seqcol

    # --- The web UI bundle --------------------------------------------------
    web_dist = resolve_web_dist(web_dist)
    build_info = read_build_info(web_dist)
    capabilities = _capabilities(mode)
    # /ping reads these off app.state.
    app.state.capabilities = capabilities

    from refgenie.server.routers.ping import router as ping_router

    app.include_router(ping_router)

    @app.get("/service-info", tags=["Default"])
    async def service_info():
        """GA4GH-style discovery document, and the web UI's bootstrap."""
        services = () if is_local else seqcol.services()
        return _service_info_doc(mode, root_path, build_info, capabilities, services)

    if not is_local:
        _mount_mcp(app, rg)

    if not is_local:
        # The llmstxt.org index for AI agents; points at /mcp.
        llms_txt = (Path(__file__).parent / "llms.txt").read_text()

        @app.get("/llms.txt", include_in_schema=False)
        def get_llms_txt() -> PlainTextResponse:
            return PlainTextResponse(llms_txt, media_type="text/plain; charset=utf-8")

    # Last: the SPA catch-all.
    mount_spa(app, web_dist, shell=is_local)

    return app


def create_local_app() -> FastAPI:
    """Factory entry point: uvicorn refgenie.server.main:create_local_app --factory"""
    return create_app(mode=APP_MODE_LOCAL)


def __getattr__(name: str):
    """Build the module-level server app on first access.

    Lazy because server mode needs the ``server`` extra, which ``refgenie dash``
    installs lack.
    """
    if name == "app":
        application = create_app()
        globals()["app"] = application
        return application
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def run_server(port: int, reload: bool = False):
    """Serve the public API on all interfaces."""
    uvicorn.run("refgenie.server.main:app", host="0.0.0.0", port=port, reload=reload)


def run_local(port: int):
    """Serve the local management app on the loopback interface only.

    Passing an app object pins uvicorn to one worker, which the in-memory
    JobManager requires.
    """
    uvicorn.run(
        create_app(mode=APP_MODE_LOCAL),
        host="127.0.0.1",
        port=port,
        log_level="warning",
        access_log=False,
    )
