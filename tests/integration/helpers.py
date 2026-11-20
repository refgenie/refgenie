"""Plain helpers for the integration tier, shared by its conftest and test modules.

Fixtures stay in ``conftest.py``; anything a test module imports lives here.
"""

import functools
import http.server
import os
import subprocess
import threading
import time
from pathlib import Path
from urllib.parse import urlparse

from tests.support.factories import TESTS_DATA_DIR, register_fasta
from tests.support.processes import run_cli

# Configuration from environment (set by test-integration.sh) or defaults
TEST_DB_URL = os.getenv(
    "TEST_DB_URL", "postgresql+psycopg://testuser:testpass@localhost:5433/refgenie_test"
)


def parse_db_url(url: str) -> dict:
    """Parse a SQLAlchemy database URL into config components."""
    parsed = urlparse(url.replace("postgresql+psycopg://", "postgresql://"))
    return {
        "type": "postgresql",
        "host": parsed.hostname or "localhost",
        "port": parsed.port or 5432,
        "user": parsed.username or "testuser",
        "password": parsed.password or "testpass",
        "name": parsed.path.lstrip("/") or "refgenie_test",
    }


def write_db_config(config_path: Path, db_url: str) -> None:
    """Write a database config YAML file for CLI usage."""
    db = parse_db_url(db_url)
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(
        f"type: {db['type']}\n"
        f"host: {db['host']}\n"
        f"port: {db['port']}\n"
        f"user: {db['user']}\n"
        f"password: {db['password']}\n"
        f"name: {db['name']}\n"
    )


def run_refgenie(
    *args: str,
    env: dict | None = None,
    check: bool = True,
) -> subprocess.CompletedProcess:
    """Run refgenie CLI command via subprocess.

    All commands are routed through the pydantic-settings CLI entry point.

    Args:
        args: CLI arguments (e.g., "init", "-f", "/path/to/genomes").
        env: Environment variables to set (merged with os.environ).
        check: If True, raise CalledProcessError on non-zero exit.

    Returns:
        CompletedProcess with stdout and stderr.
    """
    return run_cli(*args, env=env, check=check)


def normalize_fasta(text: str) -> dict[str, str]:
    """Parse FASTA into {header: sequence} dict, ignoring line wrapping."""
    sequences = {}
    header = None
    seq_parts = []
    for line in text.strip().splitlines():
        if line.startswith(">"):
            if header is not None:
                sequences[header] = "".join(seq_parts)
            header = line
            seq_parts = []
        else:
            seq_parts.append(line.strip())
    if header is not None:
        sequences[header] = "".join(seq_parts)
    return sequences


def create_client_env(tmp_path: Path, *server_urls: str) -> dict:
    """Create a fresh client environment subscribed to one or more servers.

    Uses library API instead of CLI subprocesses for faster setup.
    """
    from refgenie import Refgenie

    client_home = tmp_path / "client"
    client_home.mkdir(exist_ok=True)

    db_path = client_home / "db"
    db_config = client_home / "refgenie_db_config.yaml"
    db_config.write_text(f"type: sqlite\npath: {db_path}\n")

    genome_folder = client_home / "genomes"
    genome_folder.mkdir(exist_ok=True)

    archive_folder = client_home / "archives"
    archive_folder.mkdir(exist_ok=True)

    env = {
        "REFGENIE_HOME_PATH": str(client_home),
        "REFGENIE_DB_CONFIG_PATH": str(db_config),
        "REFGENIE_GENOME_FOLDER": str(genome_folder),
        "REFGENIE_GENOME_STAGE_FOLDER": str(archive_folder),
        "NO_COLOR": "1",
        "TERM": "dumb",
        "COLUMNS": "500",
    }

    # Init and subscribe via library API (avoids 2 subprocess calls per invocation)
    from sqlmodel import create_engine

    engine = create_engine(f"sqlite:///{db_path}")
    rg = Refgenie(database_engine=engine, suppress_migrations=True)
    rg.database.init(genome_folder=genome_folder)
    register_fasta(rg, TESTS_DATA_DIR)
    rg.servers.subscribe(server_urls=list(server_urls))
    engine.dispose()

    return env


class SlowHTTPServer:
    """Serve a directory over HTTP, waiting ``delay`` seconds before each response.

    Stands in for a remote store on S3, where each request takes a while. It
    runs in a thread of the test process, so it may only serve another
    process: a gtars call here could hold this process's GIL and stall it.
    ``first_request`` is set when the first request arrives.
    """

    def __init__(self, directory: Path, delay: float):
        self.first_request = threading.Event()
        server = self

        class Handler(http.server.SimpleHTTPRequestHandler):
            def do_GET(self):
                server.first_request.set()
                time.sleep(delay)
                super().do_GET()

            def log_message(self, *args):
                pass

        handler = functools.partial(Handler, directory=str(directory))
        self._httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
        self.url = f"http://127.0.0.1:{self._httpd.server_port}"

    def __enter__(self) -> "SlowHTTPServer":
        threading.Thread(target=self._httpd.serve_forever, daemon=True).start()
        return self

    def __exit__(self, *exc) -> None:
        self._httpd.shutdown()
        self._httpd.server_close()


def write_store_catalog(dest: Path, stores: dict[str, str]) -> None:
    """Write a publish catalog at ``dest`` that lists ``stores`` (name -> URL) and nothing else.

    The catalog is built beside ``dest`` and moved into place in one step, so
    a server importing ``dest`` on a schedule never reads half a file.
    """
    from sqlmodel import create_engine

    from refgenie import Refgenie
    from refgenie.db.tables import StoreType
    from refgenie.publish_catalog import export_publish_catalog

    build = dest.parent / f"{dest.stem}-build"
    build.mkdir(parents=True, exist_ok=True)
    db_path = build / "db.sqlite"
    db_path.unlink(missing_ok=True)
    engine = create_engine(f"sqlite:///{db_path}")
    rg = Refgenie(database_engine=engine, suppress_migrations=True)
    rg.database.init(genome_folder=build / "genomes", genome_stage_folder=build / "stage")
    for priority, (name, url) in enumerate(stores.items(), start=10):
        rg.store.add(name, url, StoreType.remote, priority=priority)
    staged = build / "catalog.sqlite"
    export_publish_catalog(rg, staged, "https://example.org/assets")
    engine.dispose()
    os.replace(staged, dest)
