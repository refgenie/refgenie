"""
Tests for the database manager (``rgc.database``): init folder creation,
idempotency, and the empty initial config, plus Refgenie constructor coercion.
"""

from pathlib import Path

from sqlmodel import select

from refgenie import Refgenie
from refgenie.db.tables import Configuration


# ---------------------------------------------------------------------------
# rgc.database.init(): folder creation, idempotency, and empty config (unit tier)
#
# Covers the regression where REFGENIE_GENOME_STAGE_FOLDER was not auto-created
# at init.
# ---------------------------------------------------------------------------


def test_init_creates_genome_and_stage_folders(engine, tmp_path):
    genome_folder = tmp_path / "g"
    stage_folder = tmp_path / "s"
    assert not genome_folder.exists()
    assert not stage_folder.exists()

    r = Refgenie(database_engine=engine, suppress_migrations=True)
    r.database.init(genome_folder=genome_folder, genome_stage_folder=stage_folder)

    assert genome_folder.is_dir()
    assert stage_folder.is_dir()


def test_init_idempotent_when_folders_exist(engine, tmp_path):
    genome_folder = tmp_path / "g"
    stage_folder = tmp_path / "s"
    genome_folder.mkdir(parents=True)
    stage_folder.mkdir(parents=True)

    r = Refgenie(database_engine=engine, suppress_migrations=True)
    # Should not raise
    r.database.init(genome_folder=genome_folder, genome_stage_folder=stage_folder)

    assert genome_folder.is_dir()
    assert stage_folder.is_dir()


def test_init_no_stage_folder_when_none(engine, tmp_path, monkeypatch):
    # Force config.genome_stage_folder to None so init() truly receives no stage folder.
    from refgenie import config as config_module

    monkeypatch.setattr(config_module.config, "genome_stage_folder", None)

    genome_folder = tmp_path / "g"
    stage_folder = tmp_path / "s"
    assert not genome_folder.exists()

    r = Refgenie(database_engine=engine, suppress_migrations=True)
    r.database.init(genome_folder=genome_folder, genome_stage_folder=None)

    assert genome_folder.is_dir()
    # With no stage folder configured, init() must not create one.
    assert not stage_folder.exists()


def test_database_init_idempotent(engine, tmp_path):
    """Regression: rgc.database.init is idempotent.

    The persistent build catalog (nightly registry builds) re-runs init on every
    run. Configuration.version is unique, so a naive re-insert raised IntegrityError
    -- logged as an alarming ERROR that also masked genuine init failures.
    init now skips the insert when a Configuration row already exists.
    """
    r = Refgenie(database_engine=engine, suppress_migrations=True)
    genome_folder = tmp_path / "genomes"
    genome_folder.mkdir()

    # First init creates the Configuration row.
    assert r.database.init(genome_folder=genome_folder) is True
    # Re-init on an already-initialized backend is a no-op success, not an
    # IntegrityError, and does not insert a duplicate row.
    assert r.database.init(genome_folder=genome_folder) is True

    with r.database._database_session as session:
        assert len(session.exec(select(Configuration)).all()) == 1


def test_initialization_produces_empty_config(engine, tmp_path):
    """After init, no asset classes or recipes are registered."""
    r = Refgenie(database_engine=engine, suppress_migrations=True)
    r.database.init(genome_folder=tmp_path / "genomes")
    assert len(r.recipe.list_all()) == 0
    assert len(r.asset_class.list_all()) == 0


def test_bare_init_never_writes_under_home(engine, tmp_path):
    """A bare init() must not resolve into the user's real ~/.refgenie.

    ``needs_migration`` can call a bare ``init()`` on its own, so the
    default has to be safe rather than fatal. The conftest redirect plus the
    autouse ``_isolated_default_genome_folder`` fixture are what make it so.
    """
    r = Refgenie(database_engine=engine, suppress_migrations=True)
    r.database.init()
    assert Path.home() / ".refgenie" not in r.genome_folder.parents
    assert r.genome_folder != Path.home() / ".refgenie" / "genomes"


class TestRefgenieInit:
    """Refgenie constructor input coercion."""

    def test_accepts_str_database_config_path(self, tmp_path):
        """Refgenie(database_config_path=<str>) coerces to Path without AttributeError."""
        config_path = tmp_path / "refgenie_db_config.yaml"
        assert Refgenie(database_config_path=str(config_path), suppress_migrations=True) is not None
        assert Refgenie(database_config_path=config_path, suppress_migrations=True) is not None
