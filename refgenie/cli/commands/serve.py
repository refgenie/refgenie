"""The `serve` and `dash` commands: models and handlers."""

import os
from typing import Literal

from pydantic import AliasChoices, BaseModel, Field

from refgenie.cli.commands.helpers import note


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

    ``refgenie.server.main`` guards only the *dash* extra's imports, so
    ``serve`` checks the ``server`` extra (apscheduler, mcp) here.
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

    note(
        f"Starting refgenie server on port {cmd.port}. "
        f"Reload mode: {'enabled' if cmd.reload else 'disabled'}."
    )
    run_server(port=cmd.port, reload=cmd.reload)


def handle_dash(cmd, refgenie) -> None:
    _require_extra("dash", "fastapi", "uvicorn")

    # LocalSecuritySettings reads this when run_local builds the app.
    if cmd.bridge is not None:
        os.environ["REFGENIE_BRIDGE_MODE"] = cmd.bridge

    from refgenie.server.main import run_local

    # 127.0.0.1, not localhost: the app binds loopback and guards the Host header.
    # Always one worker: job state lives in this process's memory.
    url = f"http://127.0.0.1:{cmd.port}"
    _open_browser_quietly(url)
    note(f"refgenie web UI at {url} (Ctrl+C to stop)")
    run_local(port=cmd.port)


def _open_browser_quietly(url: str) -> None:
    """Open ``url`` in a new tab without the browser writing to our terminal.

    fds 1 and 2 point at /dev/null while the browser spawns, so its startup
    noise stays out of the dash output.
    """
    from webbrowser import open_new_tab

    saved = [os.dup(1), os.dup(2)]
    devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull, 1)
        os.dup2(devnull, 2)
        open_new_tab(url)
    finally:
        os.dup2(saved[0], 1)
        os.dup2(saved[1], 2)
        for fd in (*saved, devnull):
            os.close(fd)
