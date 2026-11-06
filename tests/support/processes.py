"""Launch the CLI as a subprocess, pick a free port and wait for a server."""

import os
import socket
import subprocess
import sys
import time

#: The one way to launch the Python CLI in a subprocess.
PY_CLI_ARGV = [
    sys.executable,
    "-c",
    "from refgenie.cli.main import main; import sys; main(sys.argv[1:])",
]


def cli_argv(*args: str) -> list[str]:
    """``PY_CLI_ARGV`` plus ``args``, as a fresh list."""
    return PY_CLI_ARGV + [str(a) for a in args]


def _merged_env(env: dict | None) -> dict:
    """``os.environ`` updated with ``env``; a ``None`` value *unsets* the name."""
    full_env = os.environ.copy()
    for key, value in (env or {}).items():
        if value is None:
            full_env.pop(key, None)
        else:
            full_env[key] = value
    return full_env


def run_cli(
    *args: str,
    env: dict | None = None,
    check: bool = True,
    timeout: int = 120,
    text: bool = True,
) -> subprocess.CompletedProcess:
    """Run the refgenie CLI in a subprocess, capturing output.

    ``env`` is merged into ``os.environ``. ``check`` raises
    ``CalledProcessError`` *after* the output is captured (rather than using
    subprocess's own ``check=``), because that is what ``run_refgenie``'s
    callers depend on.
    """
    cmd = cli_argv(*args)
    result = subprocess.run(
        cmd, capture_output=True, text=text, env=_merged_env(env), timeout=timeout
    )
    if check and result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd, result.stdout, result.stderr)
    return result


def popen_cli(*args: str, env: dict | None = None, text: bool = False) -> subprocess.Popen:
    """Launch the refgenie CLI as a background subprocess with piped output.

    ``env`` is merged into ``os.environ``.
    """
    return subprocess.Popen(
        cli_argv(*args),
        env=_merged_env(env),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=text,
    )


def find_free_port() -> int:
    """Find a free TCP port on localhost."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def wait_for_server(host: str, port: int, timeout: float = 30.0) -> None:
    """Block until ``host:port`` accepts a connection, or raise RuntimeError."""
    start = time.time()
    while time.time() - start < timeout:
        try:
            with socket.create_connection((host, port), timeout=0.5):
                return
        except (ConnectionRefusedError, OSError):
            time.sleep(0.2)
    raise RuntimeError(f"Server failed to start at {host}:{port}")
