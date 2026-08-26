# Single-database configuration to manage refgenie database migrations

This directory contains [alembic](https://alembic.sqlalchemy.org/en/latest/) migrations for the refgenie database. The migrations are used to manage the database schema changes over time.

`1a0e5c7b9d21_initial_schema.py` is the baseline: it creates the entire current schema from `SQLModel.metadata`. It replaced an earlier chain of five revisions that no longer applied from base. There is no upgrade path from those revisions — databases predating the baseline are rebuilt, not migrated.

## How to create a new migration

Autogenerate against a database that is already at `head` (the migration is the diff between that database and `SQLModel.metadata`, so a stale database produces a wrong diff):

```bash
task alembic-revision MESSAGE="Add new table" DB_CONN_STR="postgresql://postgres:mysecretpassword@localhost:5432/postgres"
```

This writes a new file to `refgenie/db/migrations/versions/`.

### Important follow-up steps

1. Inspect the generated file. Autogenerate is a starting point, not an answer — check that `upgrade()` and `downgrade()` are inverses and that nothing was emitted for objects alembic cannot see.
2. **Update `TARGET_ALEMBIC_VERSION` in `refgenie/const.py` to the new revision id.** This is the single most important step. `Refgenie._create_db_and_tables` stamps new databases with that constant and `migrate_db()` upgrades to it, so if it does not name the chain head, every revision after it silently never runs. This has already gone wrong once.
3. Verify the result rather than assuming it:

   ```bash
   # applies from an empty database, and reverses
   alembic -x database_connection_string="sqlite:////tmp/check.db" upgrade head
   alembic -x database_connection_string="sqlite:////tmp/check.db" downgrade base
   ```

   Then confirm the migrated schema matches the models, using
   `alembic.autogenerate.compare_metadata` against `SQLModel.metadata` — it
   should report no differences. Do not eyeball it.
