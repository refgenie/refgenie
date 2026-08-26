"""
Integration test fixtures for Refgenie API.

Uses FastAPI TestClient with PostgreSQL test database.
All services are started/stopped by ./tests/scripts/test-integration.sh

Prerequisites:
    ./tests/scripts/test-integration.sh

Or manually:
    1. Start services: ./tests/scripts/services.sh start
    2. Run tests: RUN_INTEGRATION_TESTS=true pytest tests/integration/
    3. Stop services: ./tests/scripts/services.sh stop
"""

import os
import subprocess
import sys
from pathlib import Path
from urllib.parse import urlparse

import pytest

from tests.helpers import (
    DEMO_FASTA,
    RCRSD_FASTA,
    TESTS_DATA_DIR,
    find_free_port,
    make_server_client,
    popen_cli,
    register_fasta,
    requires_server,
    run_cli,
    wait_for_server,
)

GENOME_ALIAS = "rCRSd"

# Skip all integration tests unless explicitly enabled
pytestmark = pytest.mark.skipif(
    os.getenv("RUN_INTEGRATION_TESTS") != "true",
    reason="Integration tests disabled. Run ./tests/scripts/test-integration.sh",
)

# Configuration from environment (set by test-integration.sh) or defaults
TEST_DB_URL = os.getenv(
    "TEST_DB_URL", "postgresql+psycopg://testuser:testpass@localhost:5433/refgenie_test"
)
DATA_CHANNEL_URL = os.getenv("DATA_CHANNEL_URL", "http://localhost:18765")


# =============================================================================
# CLI Subprocess Helpers
# =============================================================================


def _parse_db_url(url: str) -> dict:
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


def _write_db_config(config_path: Path, db_url: str) -> None:
    """Write a database config YAML file for CLI usage."""
    db = _parse_db_url(db_url)
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


@pytest.fixture(scope="session")
def cli_home(tmp_path_factory):
    """Session-scoped temp directory serving as REFGENIE_HOME_PATH for CLI tests."""
    return tmp_path_factory.mktemp("cli_home")


@pytest.fixture(scope="session")
def cli_env(cli_home):
    """Environment variables for running refgenie CLI commands.

    Sets up a temp REFGENIE_HOME_PATH with a database config YAML pointing
    to the test PostgreSQL instance, so the CLI subprocess finds the test DB.
    """
    db_config_path = cli_home / "refgenie_db_config.yaml"
    _write_db_config(db_config_path, TEST_DB_URL)

    genome_folder = cli_home / "genomes"
    genome_folder.mkdir(exist_ok=True)

    archive_folder = cli_home / "archives"
    archive_folder.mkdir(exist_ok=True)

    return {
        "REFGENIE_HOME_PATH": str(cli_home),
        "REFGENIE_DB_CONFIG_PATH": str(db_config_path),
        "REFGENIE_GENOME_FOLDER": str(genome_folder),
        "REFGENIE_GENOME_STAGE_FOLDER": str(archive_folder),
        "RUN_INTEGRATION_TESTS": "true",
        # Disable Rich formatting in subprocess output
        "NO_COLOR": "1",
        "TERM": "dumb",
        "COLUMNS": "500",
    }


@pytest.fixture(scope="session")
def cli_build_config(cli_env, tmp_path_factory):
    """Build a real genome asset via library API for session fixtures.

    Uses the library API (not CLI subprocesses) for setup speed; CLI tests
    still exercise the CLI in their actual test assertions.
    """
    from sqlmodel import create_engine
    from refgenie import Refgenie

    genome_folder = cli_env["REFGENIE_GENOME_FOLDER"]
    stage_folder = cli_env["REFGENIE_GENOME_STAGE_FOLDER"]

    engine = create_engine(TEST_DB_URL)
    rg = Refgenie(database_engine=engine, suppress_migrations=True)
    rg.init(genome_folder=Path(genome_folder), genome_stage_folder=Path(stage_folder))
    register_fasta(rg, TESTS_DATA_DIR)

    rg.genome.initialize_genome(
        fasta_file_path=RCRSD_FASTA,
        alias_names=[GENOME_ALIAS],
        description="Simple test genome",
    )
    rg.build_asset(
        recipe_name="fasta",
        genome_name=GENOME_ALIAS,
        asset_group_name="fasta",
        asset_name="default",
    )

    genome_digest = rg.alias.resolve(GENOME_ALIAS)
    engine.dispose()

    return {
        "env": cli_env,
        "genome_folder": genome_folder,
        "genome_digest": genome_digest,
        "genome_name": GENOME_ALIAS,
    }


@pytest.fixture
def test_engine():
    """Function-scoped engine for test isolation."""
    from sqlmodel import create_engine

    engine = create_engine(TEST_DB_URL)
    yield engine
    engine.dispose()


@pytest.fixture(autouse=True)
def clean_database(test_engine, request):
    """
    Reset database to clean state before each test.

    Uses the 'shared_state' marker to skip cleaning for tests that share
    fixtures. Apply @pytest.mark.shared_state to test classes that use
    session-scoped or class-scoped fixtures.
    """
    if request.node.get_closest_marker("shared_state"):
        yield
        return

    from sqlmodel import SQLModel

    # Drop all tables
    SQLModel.metadata.drop_all(test_engine)
    # Recreate tables
    SQLModel.metadata.create_all(test_engine)

    yield


@pytest.fixture
def refgenie_instance(test_engine, tmp_path):
    """Fresh refgenie instance for each test."""
    from refgenie import Refgenie

    rg = Refgenie(database_engine=test_engine, suppress_migrations=True)
    rg.init(genome_folder=tmp_path / "genomes")
    return rg


def build_rcrsd_asset(rg, fixtures_path: Path) -> None:
    """Initialize genome and build the standard rCRSd test asset."""
    register_fasta(rg, fixtures_path)
    rg.genome.initialize_genome(
        fasta_file_path=fixtures_path / "rCRSd.fa",
        alias_names=["rCRSd"],
        description="rCRSd genome",
    )
    rg.build_asset(
        recipe_name="fasta",
        genome_name="rCRSd",
        asset_group_name="fasta",
        asset_name="test",
    )


@pytest.fixture
def client(test_engine, tmp_path):
    """TestClient over the REAL refgenie server app (create_app) with the
    test database. No hand-rolled routes -- route/schema drift in
    refgenie/server/ shows up here."""
    requires_server()
    from refgenie import Refgenie

    rg = Refgenie(database_engine=test_engine, suppress_migrations=True)
    rg.init(genome_folder=tmp_path / "genomes")
    with make_server_client(rg) as c:
        yield c


# =============================================================================
# Multi-genome and serve fixtures for bulk pull tests
# =============================================================================

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


def start_refgenie_serve(env: dict, port: int) -> subprocess.Popen:
    """Start `refgenie serve` as a subprocess on the given port.

    Uses the real server code (full OpenAPI spec), so clients can be pointed at
    it exactly as they would a production server. Raises RuntimeError with the
    captured subprocess output if the server never comes up.
    """
    proc = popen_cli("serve", "-p", str(port), env=env)
    try:
        wait_for_server("127.0.0.1", port, timeout=30)
    except RuntimeError:
        proc.terminate()
        stdout, stderr = proc.communicate(timeout=5)
        raise RuntimeError(
            f"Server failed to start at 127.0.0.1:{port}\n"
            f"STDOUT: {stdout.decode()[-500:]}\n"
            f"STDERR: {stderr.decode()[-500:]}"
        )
    return proc


def stop_refgenie_serve(proc: subprocess.Popen) -> None:
    """Terminate a `refgenie serve` subprocess, killing it if it hangs."""
    proc.terminate()
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait()


def build_server_env(tmp_dir: Path, genomes: list[dict]) -> dict:
    """Build a complete server environment with multiple genomes + fasta assets.

    Each genome spec is a dict with keys:
        - fasta: Path to the FASTA file
        - alias: str genome alias
        - override: list[str] | None -- serving_modes_override

    Creates a single SQLite DB, one Refgenie instance, registers fasta once,
    then loops over each genome spec to init, build, optionally override, and stage.

    Returns dict with "engine" (for teardown), "env" (CLI/subprocess environment
    variables) and a "genomes" dict keyed by alias containing genome_digest and
    asset_digest for each genome.
    """
    from sqlmodel import Session, SQLModel, create_engine

    from refgenie import Refgenie
    from refgenie.db.tables import Asset

    db_path = tmp_dir / "refgenie.db"
    genome_folder = tmp_dir / "genomes"
    stage_folder = tmp_dir / "archives"
    genome_folder.mkdir(exist_ok=True)
    stage_folder.mkdir(exist_ok=True)

    engine = create_engine(f"sqlite:///{db_path}")
    SQLModel.metadata.create_all(engine)

    rg = Refgenie(database_engine=engine, suppress_migrations=True)
    rg.init(genome_folder=genome_folder, genome_stage_folder=stage_folder)
    register_fasta(rg, TESTS_DATA_DIR)

    genomes_info = {}
    for spec in genomes:
        fasta_path = spec["fasta"]
        alias = spec["alias"]
        override = spec.get("override")

        rg.genome.initialize_genome(
            fasta_file_path=fasta_path,
            alias_names=[alias],
            description=f"{alias} test genome",
        )
        rg.build_asset(
            recipe_name="fasta",
            genome_name=alias,
            asset_group_name="fasta",
            asset_name="default",
        )

        genome_digest = rg.alias.resolve(alias)
        asset = rg.asset.get(
            genome_name=alias,
            asset_group_name="fasta",
            asset_name="default",
        )
        asset_digest = asset.digest

        if override:
            from sqlalchemy import update as sa_update

            with Session(engine) as session:
                session.exec(
                    sa_update(Asset)
                    .where(Asset.digest == asset_digest)
                    .values(serving_modes_override=override)
                )
                session.commit()
            asset.serving_modes_override = override

        rg.stage.create(
            asset=asset,
            genome_folder=genome_folder,
            genome_stage_folder=stage_folder,
        )

        genomes_info[alias] = {
            "genome_digest": genome_digest,
            "asset_digest": asset_digest,
        }

    # Create env dict for CLI subprocess usage
    db_config = tmp_dir / "refgenie_db_config.yaml"
    db_config.write_text(f"type: sqlite\npath: {db_path}\n")

    env = {
        "REFGENIE_HOME_PATH": str(tmp_dir),
        "REFGENIE_DB_CONFIG_PATH": str(db_config),
        "REFGENIE_GENOME_FOLDER": str(genome_folder),
        "REFGENIE_GENOME_STAGE_FOLDER": str(stage_folder),
        "NO_COLOR": "1",
        "TERM": "dumb",
        "COLUMNS": "500",
    }

    return {
        "engine": engine,
        "genomes": genomes_info,
        "env": env,
    }


@pytest.fixture(scope="session")
def multi_remote_servers(tmp_path_factory):
    """Two real `refgenie serve` subprocesses: A has rCRSd, B has demo.

    Both assets use the default fasta asset class serving modes
    (``[file, archive]``), so server B doubles as the both-modes fixture.
    """
    tmp_a = tmp_path_factory.mktemp("server_a")
    info_a = build_server_env(tmp_a, [{"fasta": RCRSD_FASTA, "alias": "rCRSd"}])

    tmp_b = tmp_path_factory.mktemp("server_b")
    info_b = build_server_env(tmp_b, [{"fasta": DEMO_FASTA, "alias": "demo"}])

    port_a = find_free_port()
    port_b = find_free_port()

    proc_a = start_refgenie_serve(info_a["env"], port_a)
    proc_b = start_refgenie_serve(info_b["env"], port_b)

    yield {
        "server_a": {**info_a, "url": f"http://127.0.0.1:{port_a}"},
        "server_b": {**info_b, "url": f"http://127.0.0.1:{port_b}"},
    }

    stop_refgenie_serve(proc_a)
    stop_refgenie_serve(proc_b)
    info_a["engine"].dispose()
    info_b["engine"].dispose()


@pytest.fixture(scope="session")
def file_mode_pull_server(tmp_path_factory):
    """Real `refgenie serve` subprocess with a file-only (override) fasta asset."""
    tmp = tmp_path_factory.mktemp("file_pull")
    info = build_server_env(
        tmp,
        [{"fasta": RCRSD_FASTA, "alias": "rCRSd", "override": ["file"]}],
    )

    port = find_free_port()
    proc = start_refgenie_serve(info["env"], port)

    yield {**info, "url": f"http://127.0.0.1:{port}"}

    stop_refgenie_serve(proc)
    info["engine"].dispose()


def _create_client_env(tmp_path: Path, *server_urls: str) -> dict:
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
    rg.init(genome_folder=genome_folder)
    register_fasta(rg, TESTS_DATA_DIR)
    rg.configuration.subscribe(server_urls=list(server_urls))
    engine.dispose()

    return env


@pytest.fixture(scope="session")
def cli_multi_genome_config(cli_build_config):
    """Build a second genome (demo) for multi-genome tests via library API.

    Depends on cli_build_config which already built rCRSd.
    """
    from sqlmodel import create_engine
    from refgenie import Refgenie

    engine = create_engine(TEST_DB_URL)
    rg = Refgenie(database_engine=engine, suppress_migrations=True)
    rg.init(genome_folder=Path(cli_build_config["genome_folder"]))
    register_fasta(rg, TESTS_DATA_DIR)

    try:
        rg.genome.initialize_genome(
            fasta_file_path=DEMO_FASTA,
            alias_names=["demo"],
            description="demo genome",
        )
    except ValueError:
        # Already exists from a previous run
        pass

    try:
        rg.build_asset(
            recipe_name="fasta",
            genome_name="demo",
            asset_group_name="fasta",
            asset_name="default",
        )
    except Exception:
        pass  # Asset may already exist

    engine.dispose()

    return cli_build_config


@pytest.fixture(scope="session")
def cli_archived_config(cli_build_config):
    """Stage the built asset via library API.

    Session-scoped fixture that creates archives.
    """
    from sqlmodel import create_engine
    from refgenie import Refgenie

    env = cli_build_config["env"]
    archive_folder = env["REFGENIE_GENOME_STAGE_FOLDER"]
    genome_folder = Path(cli_build_config["genome_folder"])

    engine = create_engine(TEST_DB_URL)
    rg = Refgenie(database_engine=engine, suppress_migrations=True)
    # Configuration already exists with genome_stage_folder from cli_build_config

    asset = rg.asset.get(
        genome_name=GENOME_ALIAS,
        asset_group_name="fasta",
        asset_name="default",
    )
    rg.stage.create(
        asset=asset,
        genome_folder=genome_folder,
        genome_stage_folder=Path(archive_folder),
    )

    # Get archive info
    from sqlmodel import Session, select
    from sqlalchemy.orm import selectinload
    from refgenie.db.tables import StagedAsset, Asset as AssetModel, Configuration

    with Session(engine) as session:
        config = session.exec(select(Configuration)).first()
        staged = session.exec(
            select(StagedAsset)
            .options(selectinload(StagedAsset.asset).selectinload(AssetModel.asset_group))
            .where(StagedAsset.mode == "archive")
        ).first()
        archive_digest = staged.asset_digest if staged else None
        if staged and config and config.genome_stage_folder:
            from refgenie.utils.staging import staged_archive_path

            archive_path = str(
                staged_archive_path(
                    config.genome_stage_folder,
                    staged.asset.asset_group.genome_digest,
                    staged.asset.asset_group.name,
                    staged.asset.digest,
                )
            )
        else:
            archive_path = None
        archive_size = staged.tarball_size if staged else None

    engine.dispose()

    return {
        **cli_build_config,
        "archive_digest": archive_digest,
        "archive_path": archive_path,
        "archive_size": archive_size,
        "archive_folder": archive_folder,
    }


@pytest.fixture(scope="session")
def cli_multi_archived_config(cli_archived_config, cli_multi_genome_config):
    """Archive both genomes via library API.

    cli_archived_config already archived rCRSd.
    cli_multi_genome_config already built demo.
    This archives demo too.
    """
    from sqlmodel import create_engine, Session, select
    from refgenie import Refgenie
    from refgenie.db.tables import StagedAsset

    genome_folder = Path(cli_archived_config["genome_folder"])
    archive_folder = Path(cli_archived_config["archive_folder"])

    engine = create_engine(TEST_DB_URL)
    rg = Refgenie(database_engine=engine, suppress_migrations=True)
    # Configuration already exists with genome_stage_folder from cli_archived_config

    try:
        asset = rg.asset.get(
            genome_name="demo",
            asset_group_name="fasta",
            asset_name="default",
        )
        rg.stage.create(
            asset=asset,
            genome_folder=genome_folder,
            genome_stage_folder=archive_folder,
        )
    except Exception:
        pass  # Staged asset may already exist

    # Verify both archives exist
    with Session(engine) as session:
        archive_count = len(
            session.exec(select(StagedAsset).where(StagedAsset.mode == "archive")).all()
        )
    engine.dispose()

    assert archive_count >= 2, (
        f"Expected at least 2 staged archive records (rCRSd + demo), found {archive_count}"
    )

    return cli_archived_config


@pytest.fixture(scope="session")
def refgenie_serve_subprocess(cli_multi_archived_config):
    """Start refgenie serve as a real subprocess.

    This tests the actual server code, not a mock.
    Both rCRSd and demo are built and archived.
    """
    env = cli_multi_archived_config["env"]
    port = find_free_port()
    server_url = f"http://127.0.0.1:{port}"

    proc = start_refgenie_serve(env, port)

    yield {
        **cli_multi_archived_config,
        "url": server_url,
        "port": port,
        "process": proc,
    }

    stop_refgenie_serve(proc)


# =============================================================================
# Store-backed serving fixtures
# =============================================================================
# These exercise refgenie's production server path: a remote RefgetStore served
# over HTTP, with the refgenie server running in store-backed mode
# (REFGENIE_REFGET_STORE_URL set). DB-backed mode is the default and is covered
# elsewhere; this keeps that coverage intact while adding the store-backed path.
#
# NOTE: gtars' RefgetStore.open_remote (Rust/PyO3) holds the GIL during HTTP
# requests, which deadlocks a Python-thread-based HTTP server. The store MUST be
# served from a separate process. We therefore serve it with `python -m
# http.server` (subprocess) and run the refgenie server as a subprocess too --
# matching the existing `refgenie serve` subprocess pattern.

# Small demo collections from the refget package's test FASTA files.
# Resolved relative to this repo's sibling `refget` source checkout.
_REFGET_TEST_FASTA_DIR = Path(__file__).resolve().parents[3] / "refget" / "test_fasta"
STORE_DEMO_FASTAS = ["base.fa", "subset.fa", "different_names.fa"]


def _serve_directory_http(directory: Path) -> tuple[subprocess.Popen, str]:
    """Serve a directory over HTTP via a `python -m http.server` subprocess.

    Returns (process, base_url). A subprocess (not a thread) is required because
    gtars' open_remote holds the GIL during HTTP requests and would deadlock a
    Python-thread HTTP server.
    """
    port = find_free_port()
    proc = subprocess.Popen(
        [sys.executable, "-m", "http.server", str(port), "--directory", str(directory)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    try:
        wait_for_server("127.0.0.1", port, timeout=15)
    except RuntimeError:
        proc.terminate()
        raise
    return proc, f"http://127.0.0.1:{port}"


def _build_demo_store(store_dir: Path) -> list[str]:
    """Build an on-disk RefgetStore from the small refget test FASTAs.

    Returns the list of collection digests (FASTA-add order).
    """
    from refget.store import RefgetStore

    if not _REFGET_TEST_FASTA_DIR.is_dir():
        raise RuntimeError(
            f"refget test_fasta dir not found at {_REFGET_TEST_FASTA_DIR}; "
            "store-backed integration tests require the refget source checkout."
        )

    store = RefgetStore.on_disk(str(store_dir))
    for fa in STORE_DEMO_FASTAS:
        store.add_sequence_collection_from_fasta(str(_REFGET_TEST_FASTA_DIR / fa))
    return [c.digest for c in store.list_collections()["results"]]


@pytest.fixture(scope="session")
def store_backed_server(tmp_path_factory):
    """Run a refgenie server in store-backed mode against an HTTP-served store.

    Builds a small on-disk RefgetStore (3 collections from refget test FASTAs),
    serves it over HTTP, then starts `refgenie serve` with
    REFGENIE_REFGET_STORE_URL pointing at the HTTP store (production path).

    Yields a dict with:
        - url: refgenie server base URL
        - store_url: HTTP URL of the served RefgetStore
        - digests: list of collection digests in the store
    """
    tmp = tmp_path_factory.mktemp("store_backed")

    # 1. Build the on-disk store and serve it over HTTP (subprocess).
    store_dir = tmp / "store"
    store_dir.mkdir()
    digests = _build_demo_store(store_dir)
    store_proc, store_url = _serve_directory_http(store_dir)

    # 2. Start refgenie serve in store-backed mode (subprocess).
    home = tmp / "home"
    home.mkdir()
    db_config = home / "refgenie_db_config.yaml"
    db_config.write_text(f"type: sqlite\npath: {home / 'db'}\n")
    genome_folder = home / "genomes"
    genome_folder.mkdir()

    env = {
        "REFGENIE_HOME_PATH": str(home),
        "REFGENIE_DB_CONFIG_PATH": str(db_config),
        "REFGENIE_GENOME_FOLDER": str(genome_folder),
        "REFGENIE_REFGET_STORE_URL": store_url,
        "NO_COLOR": "1",
        "TERM": "dumb",
        "COLUMNS": "500",
    }

    server_port = find_free_port()
    # ``None`` unsets: no stale REFGENIE config var leaks into the subprocess.
    server_env = {**env, "REFGENIE": None}

    proc = popen_cli("serve", "-p", str(server_port), env=server_env)

    try:
        wait_for_server("127.0.0.1", server_port, timeout=40)
    except RuntimeError:
        proc.terminate()
        store_proc.terminate()
        stdout, stderr = proc.communicate(timeout=5)
        raise RuntimeError(
            f"Store-backed server failed to start at 127.0.0.1:{server_port}\n"
            f"STDOUT: {stdout.decode()[-1000:]}\n"
            f"STDERR: {stderr.decode()[-1000:]}"
        )

    yield {
        "url": f"http://127.0.0.1:{server_port}",
        "store_url": store_url,
        "digests": digests,
    }

    # Cleanup: refgenie server first, then the HTTP store server.
    proc.terminate()
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait()
    store_proc.terminate()
    try:
        store_proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        store_proc.kill()
        store_proc.wait()


@pytest.fixture(scope="session")
def bulk_pull_client_env(refgenie_serve_subprocess, tmp_path_factory):
    """Shared client environment for BulkPull tests.

    Created once per session for speed. Tests use --force to avoid conflicts
    from shared state.
    """
    tmp = tmp_path_factory.mktemp("bulk_client")
    return _create_client_env(tmp, refgenie_serve_subprocess["url"])


@pytest.fixture(scope="session")
def data_channel_url():
    """URL of the HTTP data channel server (started by test-integration.sh)."""
    return DATA_CHANNEL_URL


# =============================================================================
# Shared fixtures for read-only tests
# =============================================================================
# Session-scoped fixtures that build rCRSd ONCE for all read-only tests.
# Uses a separate SQLite database to avoid conflicting with the session
# PostgreSQL database used by CLI tests.


@pytest.fixture(scope="session")
def shared_test_engine(tmp_path_factory):
    """
    Session-scoped engine for shared read-only tests.

    Uses a separate SQLite database (not the shared PostgreSQL) to avoid
    conflicting with session CLI fixtures. This means read-only library
    tests build once per session instead of once per class.
    """
    from sqlmodel import create_engine, SQLModel

    db_dir = tmp_path_factory.mktemp("shared_test_db")
    db_path = db_dir / "shared_test.db"
    engine = create_engine(f"sqlite:///{db_path}")
    SQLModel.metadata.create_all(engine)
    yield engine
    engine.dispose()


@pytest.fixture(scope="session")
def shared_rcrsd_instance(shared_test_engine, fixtures_path, tmp_path_factory):
    """
    Session-scoped refgenie instance with pre-built rCRSd asset.

    Builds the rCRSd FASTA asset ONCE per session and shares it across all
    read-only tests. Uses a separate SQLite database to avoid conflicting
    with the session PostgreSQL fixtures.
    """
    from refgenie import Refgenie

    tmp_dir = tmp_path_factory.mktemp("shared_genomes")

    rg = Refgenie(database_engine=shared_test_engine, suppress_migrations=True)
    rg.init(genome_folder=tmp_dir / "genomes")
    register_fasta(rg, fixtures_path)

    # Initialize genome and build the rCRSd asset once
    rg.genome.initialize_genome(
        fasta_file_path=fixtures_path / "rCRSd.fa",
        alias_names=["rCRSd"],
        description="rCRSd genome",
    )
    rg.build_asset(
        recipe_name="fasta",
        genome_name="rCRSd",
        asset_group_name="fasta",
        asset_name="test",
    )

    return rg


@pytest.fixture(scope="session")
def shared_client(shared_rcrsd_instance):
    """
    Session-scoped TestClient over the REAL server app wrapping
    shared_rcrsd_instance (same engine, same genome folder). Used by
    read-only tests that need both a library instance and HTTP endpoints
    pointing to the same data -- /v4/genomes/{digest} serves the real JSON
    API over the built asset.
    """
    requires_server()
    with make_server_client(shared_rcrsd_instance) as c:
        yield c
