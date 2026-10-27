"""The `serve` and `dash` commands: models and handlers."""

import os
from typing import Literal

from pydantic import AliasChoices, BaseModel, Field

from refgenie.logger import logger


class ServeModel(BaseModel):
    """Start the production refgenie server."""

    port: int = Field(
        default=8000,
        description="Port to run the server on.",
        validation_alias=AliasChoices("p", "port"),
    )
    reload: bool = Field(
        default=False,
        description="Enable auto-reload on code changes (for development).",
        validation_alias=AliasChoices("r", "reload"),
    )


class DashModel(BaseModel):
    """Start the local refgenie web UI."""

    port: int = Field(
        default=8080,
        description="Port to run the dashboard on.",
        validation_alias=AliasChoices("p", "port"),
    )
    bridge: Literal["off", "read", "full"] | None = Field(
        default=None,
        description=(
            "Localhost-bridge mode for this run (overrides $REFGENIE_BRIDGE_MODE): "
            "'off' = no cross-origin access, 'read' = allowlisted public origins "
            "may read, 'full' = additionally allows cross-origin pull."
        ),
        validation_alias=AliasChoices("b", "bridge"),
    )


def _require_extra(extra: str, *modules: str) -> None:
    """Fail with the "install this extra" message if ``modules`` are missing.

    Both commands now go through ``refgenie.server.main``, which guards only the
    imports the *dash* extra provides. Each command still has to say which extra
    the user actually needs: ``refgenie serve`` additionally needs apscheduler
    and mcp, which live in the ``server`` extra.
    """
    from importlib import import_module

    for module in modules:
        try:
            import_module(module)
        except ImportError as e:
            raise ImportError(
                f"The '{extra}' extras are not installed. Please install refgenie "
                f"with the '{extra}' extras."
            ) from e


def handle_serve(cmd, refgenie) -> None:
    _require_extra("server", "fastapi", "uvicorn", "apscheduler", "mcp")
    from refgenie.server.main import run_server

    logger.info(
        f"Starting refgenie server on port {cmd.port}. "
        f"Reload mode: {'enabled' if cmd.reload else 'disabled'}."
    )
    run_server(port=cmd.port, reload=cmd.reload)


def handle_dash(cmd, refgenie) -> None:
    _require_extra("dash", "fastapi", "uvicorn")
    from webbrowser import open_new_tab

    # --bridge overrides the env default for this single run. Setting the env
    # var (rather than threading a parameter through run_local/create_app) works
    # because LocalSecuritySettings reads REFGENIE_BRIDGE_MODE at app
    # construction, which happens inside run_local, after this line.
    if cmd.bridge is not None:
        os.environ["REFGENIE_BRIDGE_MODE"] = cmd.bridge

    from refgenie.server.main import run_local

    # 127.0.0.1, not localhost: the local app binds the loopback address and
    # its Host-header guard is what makes that binding meaningful.
    #
    # Single worker, always. Background job state (pull/build progress) lives
    # in one process's memory, so a second worker would give the browser a
    # coin-flip 404 on the job it just submitted. Do not add a --workers flag
    # here without making jobs shared state first.
    open_new_tab(f"http://127.0.0.1:{cmd.port}")
    logger.info(
        f"Starting the refgenie web UI on port {cmd.port}. "
        "Refresh the page if it does not load automatically."
    )
    run_local(port=cmd.port)
