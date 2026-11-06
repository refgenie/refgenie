"""refgenie.db: ORM models, the alembic chain, database config and subscriptions."""

import importlib
from pathlib import Path

import pytest
from alembic.autogenerate import compare_metadata
from alembic.config import Config
from alembic.migration import MigrationContext
from alembic.script import ScriptDirectory
from sqlalchemy import create_engine, inspect
from sqlalchemy.engine.url import make_url
from sqlmodel import Session, SQLModel, select

import refgenie.config.env as config_env
import refgenie.db.tables  # noqa: F401  -- populates SQLModel.metadata
from refgenie import Refgenie
from refgenie.config.db import DatabaseType, create_default_db_config
from refgenie.const import TARGET_ALEMBIC_VERSION
from refgenie.db.migrations.utils import run_sql_migrations
from refgenie.exceptions import DatabaseMigrationError


class TestDatabaseModels:
    """refgenie.db.tables: genome/asset persistence and ORM relationships."""

    def test_genome_round_trip_by_digest(self, refgenie_session):
        """alias.resolve yields a 32-char seqcol digest; genome.get round-trips it."""
        genome_digest = refgenie_session.alias.resolve("rCRSd")
        assert genome_digest is not None
        assert len(genome_digest) == 32  # seqcol digest length

        genome = refgenie_session.genome.get(genome_digest)
        assert genome is not None
        assert genome.digest == genome_digest

    def test_asset_relationship_chain(self, refgenie_session):
        """The asset → asset_group → genome ORM chain resolves both ways."""
        genome_digest = refgenie_session.alias.resolve("rCRSd")
        asset = refgenie_session.asset.get(
            genome_digest=genome_digest, asset_group_name="fasta", asset_name="test"
        )
        assert asset is not None
        assert asset.asset_group is not None
        assert asset.asset_group.name == "fasta"
        assert asset.asset_group.genome.digest == genome_digest


class TestManagerRoot:
    """Registry root smoke coverage.

    The only place configuration.list_all, stage.list_all, and
    sources.list_channels are exercised.
    """

    def test_root_returns_lists(self, refgenie_minimal):
        r = refgenie_minimal
        assert isinstance(r.configuration.list_all(), list)
        assert isinstance(r.stage.list_all(), list)
        assert isinstance(r.sources.list_channels(), list)


# ---------------------------------------------------------------------------
# Guards that keep the alembic migration chain and the code that uses it in step.
# After adding a migration, bump TARGET_ALEMBIC_VERSION in refgenie/const.py.
# ---------------------------------------------------------------------------

VERSION_TABLE = "alembic_version"  # alembic-managed; never part of a migration's DDL
MIGRATIONS_DIR = Path(refgenie.db.tables.__file__).parent / "migrations"


def _script_directory() -> ScriptDirectory:
    config = Config()
    config.set_main_option("script_location", str(MIGRATIONS_DIR))
    return ScriptDirectory.from_config(config)


def _downgrade(db_url, revision: str) -> None:
    """``run_sql_migrations`` in reverse; there is no shipped downgrade helper."""
    from alembic import command

    config = Config()
    config.set_main_option("script_location", str(MIGRATIONS_DIR))
    config.set_main_option(
        "sqlalchemy.url", db_url.render_as_string(hide_password=False).replace("%", "%%")
    )
    config.set_main_option("path_separator", "os")
    command.downgrade(config=config, revision=revision)


class TestAlembicChain:
    """The chain head and the schema it produces must stay aligned with the models."""

    def test_target_version_is_chain_head(self):
        """
        TARGET_ALEMBIC_VERSION must name the single head of the chain. If it names
        anything else, every revision after it is dead code: new databases are
        stamped to a version they are not and needs_migration() reports
        them up to date. Nothing raises; the schema is simply wrong.
        """
        heads = _script_directory().get_heads()
        assert len(heads) == 1, (
            f"expected a single migration head, found {heads}. "
            "A branched chain means someone added a revision without rebasing onto the head."
        )
        assert TARGET_ALEMBIC_VERSION == heads[0], (
            f"TARGET_ALEMBIC_VERSION is {TARGET_ALEMBIC_VERSION!r} but the chain head is "
            f"{heads[0]!r}. Update TARGET_ALEMBIC_VERSION in refgenie/const.py to the new head, "
            "or the migrations between the two will never run."
        )

    def test_upgrade_from_base_matches_models(self, tmp_path):
        """
        `alembic upgrade head` on an empty database must reproduce
        SQLModel.metadata — tables, columns, and (via compare_metadata) types,
        nullability, indexes, unique constraints and FKs. Also catches a chain
        that cannot apply from base at all. tmp_path is a throwaway.
        """
        # Scope the comparison to refgenie's own tables. SQLModel.metadata is a
        # process-global registry: gtars' refget-store models (sequence,
        # sequencecollection, namesattr, ...) register into it too once any test
        # opens a store, so comparing against the whole metadata would report
        # those foreign tables as drift depending on test order.
        model_table_names = {
            m.__tablename__
            for m in vars(refgenie.db.tables).values()
            if isinstance(m, type)
            and issubclass(m, SQLModel)
            and getattr(m, "__table__", None) is not None
        } - {VERSION_TABLE}

        db_path = tmp_path / "migrated.db"
        run_sql_migrations(db_url=make_url(f"sqlite:///{db_path}"))

        engine = create_engine(f"sqlite:///{db_path}")
        try:
            inspector = inspect(engine)

            db_tables = set(inspector.get_table_names()) - {VERSION_TABLE}
            assert model_table_names == db_tables, (
                f"tables only in the models: {sorted(model_table_names - db_tables)}; "
                f"tables only in the migrated database: {sorted(db_tables - model_table_names)}"
            )

            for table in sorted(model_table_names):
                model_columns = {c.name for c in SQLModel.metadata.tables[table].columns}
                db_columns = {c["name"] for c in inspector.get_columns(table)}
                assert model_columns == db_columns, (
                    f"{table}: columns only in the model {sorted(model_columns - db_columns)}; "
                    f"columns only in the migrated database {sorted(db_columns - model_columns)}"
                )

            # Beyond names: types, nullability, indexes, unique constraints, FKs.
            # include_object restricts the diff to refgenie's tables so foreign
            # SQLModel tables sharing the global metadata are ignored.
            def _include_object(obj, name, type_, reflected, compare_to):
                if type_ == "table":
                    return name in model_table_names or name == VERSION_TABLE
                return True

            with engine.connect() as conn:
                context = MigrationContext.configure(
                    conn, opts={"compare_type": True, "include_object": _include_object}
                )
                diffs = [
                    d
                    for d in compare_metadata(context, SQLModel.metadata)
                    if VERSION_TABLE not in repr(d)
                ]
            assert not diffs, "migrated schema differs from SQLModel.metadata:\n" + "\n".join(
                f"  {d}" for d in diffs
            )
        finally:
            engine.dispose()

    def test_downgrade_to_base_drops_everything(self, tmp_path):
        """`alembic downgrade base` after `upgrade head` leaves no tables behind,
        and the schema can then be upgraded again."""
        url = make_url(f"sqlite:///{tmp_path / 'refgenie.db'}")
        run_sql_migrations(db_url=url)
        _downgrade(url, "base")
        engine = create_engine(url)
        try:
            assert set(inspect(engine).get_table_names()) - {VERSION_TABLE} == set()
        finally:
            engine.dispose()
        run_sql_migrations(db_url=url)

    def test_needs_migration_does_not_leak_connections(self, tmp_path):
        """needs_migration must return its pooled connection every call;
        `MigrationContext.configure(engine.connect())` once never closed it."""
        engine = create_engine(f"sqlite:///{tmp_path / 'refgenie.db'}", echo=False)
        r = Refgenie(database_engine=engine, suppress_migrations=True)
        r.database.init(genome_folder=tmp_path / "genomes")

        baseline = engine.pool.checkedout()
        for _ in range(5):
            r.database.needs_migration()

        assert engine.pool.checkedout() == baseline

    def test_failed_migration_raises(self, tmp_path):
        """A database stamped with a revision the chain does not know cannot be
        migrated; Refgenie must refuse to start rather than run on it."""
        engine = create_engine(f"sqlite:///{tmp_path / 'refgenie.db'}", echo=False)
        r = Refgenie(database_engine=engine, suppress_migrations=True)
        r.database.init(genome_folder=tmp_path / "genomes")
        with engine.begin() as conn:
            conn.exec_driver_sql("UPDATE alembic_version SET version_num = 'deadbeef0000'")
        engine.dispose()

        with pytest.raises(DatabaseMigrationError, match="deadbeef0000"):
            Refgenie(database_engine=engine)


# ---------------------------------------------------------------------------
# Configuration: database config, migration state, the legacy env-var warning,
# and server subscriptions
#
# Kept to the real behavioral contracts: the default config, the instance's
# resolved config, the migration-needed check, the import-time warning fired
# when legacy ``$REFGENIE`` is set but ``$REFGENIE_DB_CONFIG_PATH`` is not,
# and subscribe/unsubscribe (pure database operations on the Configuration
# table). (Schema/migration mechanics are covered by TestAlembicChain above,
# including the migration-check connection-leak edge.)
# ---------------------------------------------------------------------------


class TestDatabaseConfig:
    """refgenie.config.db: default and instance database configuration."""

    def test_default_db_config_is_sqlite(self):
        """The default config is SQLite with a sqlite:///…/refgenie URL."""
        config = create_default_db_config()
        assert config.type == DatabaseType.SQLITE
        assert config.url.startswith("sqlite:///")
        assert config.url.endswith("/refgenie")

    def test_instance_database_config(self, refgenie_minimal):
        """load_database_config resolves to a SQLite config with a sqlite:/// URL."""
        from refgenie.db.engine import load_database_config

        config = load_database_config()
        assert config.type == DatabaseType.SQLITE
        assert config.url.startswith("sqlite:///")


class TestMigrationCheck:
    """The catalog reports whether it needs a schema upgrade."""

    def test_freshly_initialized_catalog_needs_no_migration(self, refgenie_minimal):
        """A freshly initialized, head-stamped catalog reports no migration needed."""
        assert refgenie_minimal.database.needs_migration(log=False) is False


def test_warning_fires_when_REFGENIE_set_and_new_var_missing(monkeypatch, capsys):
    monkeypatch.setenv("REFGENIE", "/legacy/genome_config.yaml")
    monkeypatch.delenv("REFGENIE_DB_CONFIG_PATH", raising=False)
    importlib.reload(config_env)
    captured = capsys.readouterr()
    assert "REFGENIE_DB_CONFIG_PATH" in captured.err
    assert "/legacy/genome_config.yaml" in captured.err
    assert "refgenie1" not in captured.err


@pytest.mark.parametrize(
    "legacy_set, new_var_set",
    [
        pytest.param(True, True, id="new_var_set"),
        pytest.param(False, False, id="REFGENIE_unset"),
    ],
)
def test_warning_silent(monkeypatch, capsys, tmp_path, legacy_set, new_var_set):
    """No warning when the new var is set, or when legacy $REFGENIE is unset."""
    if legacy_set:
        monkeypatch.setenv("REFGENIE", "/legacy/genome_config.yaml")
    else:
        monkeypatch.delenv("REFGENIE", raising=False)
    if new_var_set:
        monkeypatch.setenv("REFGENIE_DB_CONFIG_PATH", str(tmp_path / "db.yaml"))
    else:
        monkeypatch.delenv("REFGENIE_DB_CONFIG_PATH", raising=False)
    importlib.reload(config_env)
    captured = capsys.readouterr()
    assert "WARNING" not in captured.err


S1 = "http://server1.example.com"
S2 = "http://server2.example.com"
OLD = "http://old-server.example.com"
NEW = "http://new-server.example.com"


class TestSubscribe:
    """Test server subscription management."""

    @pytest.mark.parametrize(
        "ops, expected",
        [
            ([("subscribe", S1, {})], [S1]),
            ([("subscribe", [S1, S2], {})], [S1, S2]),
            ([("subscribe", S1, {}), ("subscribe", S1, {})], [S1]),
            ([("subscribe", OLD, {}), ("subscribe", NEW, {"reset": True})], [NEW]),
            ([("subscribe", [S1, S2], {}), ("unsubscribe", [S1], {})], [S2]),
            ([("subscribe", S1, {}), ("unsubscribe", ["http://nope.example.com"], {})], [S1]),
            ([], []),
        ],
        ids=[
            "single",
            "multiple",
            "deduplicates",
            "reset",
            "unsubscribe-one",
            "unsubscribe-unknown-is-noop",
            "none-subscribed",
        ],
    )
    def test_subscription_crud(self, refgenie_minimal, ops, expected):
        """Each row applies its operations in order, then pins the full subscription set."""
        for method, arg, kwargs in ops:
            getattr(refgenie_minimal.servers, method)(arg, **kwargs)
        subs = refgenie_minimal.servers.subscriptions()
        assert sorted(subs) == sorted(expected)

    def test_subscribe_when_configuration_id_is_not_one(self, refgenie_minimal):
        """The sole Configuration row need not have id 1 (restored dump, merged db)."""
        from refgenie.db.tables import Configuration

        with Session(refgenie_minimal.database_engine) as session:
            row = session.exec(select(Configuration)).one()
            session.delete(row)
            session.commit()
            session.add(
                Configuration(
                    id=42,
                    version=row.version,
                    genome_folder=row.genome_folder,
                    genome_stage_folder=row.genome_stage_folder,
                )
            )
            session.commit()

        refgenie_minimal.servers.subscribe("http://server1.example.com")
        assert "http://server1.example.com" in refgenie_minimal.configuration.get_latest().servers
        refgenie_minimal.servers.unsubscribe(["http://server1.example.com"])
        assert refgenie_minimal.servers.subscriptions() == []
