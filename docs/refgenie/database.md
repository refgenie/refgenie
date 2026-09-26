# Database backends

Refgenie stores all asset metadata, genome records, and configuration in a database. Two backends are supported: SQLite (default) and PostgreSQL.

## SQLite (default)

SQLite requires zero configuration. When you run `refgenie init`, a SQLite database is created automatically at `~/.refgenie/refgenie`. This is the right choice for single-user workstations.

The default configuration is equivalent to:

```yaml
# ~/.refgenie/refgenie_db_config.yaml
type: sqlite
path: ~/.refgenie/refgenie
```

You do not need to create this file -- refgenie generates it automatically if it does not exist.

## PostgreSQL

For multi-user environments, server deployments, or when you need concurrent access, use PostgreSQL.

Create a database configuration file:

```yaml
# ~/.refgenie/refgenie_db_config.yaml
type: postgresql
name: refgenie
user: refgenie_user
password: your_secure_password
host: localhost
port: 5432
```

The `dialect` field is optional and defaults to `psycopg`. Refgenie uses the `psycopg` driver (PostgreSQL 3+ adapter).

## Setting the config path

Refgenie looks for the database configuration file at the path specified by the `REFGENIE_DB_CONFIG_PATH` environment variable. If this variable is not set, it defaults to `$REFGENIE_HOME_PATH/refgenie_db_config.yaml`.

```bash
export REFGENIE_DB_CONFIG_PATH=~/.refgenie/refgenie_db_config.yaml
```

## When to use which backend

| Scenario | Recommended backend |
|---|---|
| Single user on a workstation | SQLite (default) |
| Shared lab or team server | PostgreSQL |
| Production refgenie server (`refgenie serve`) | PostgreSQL |
| CI/CD or ephemeral environments | SQLite |
| Cloud deployment with multiple workers | PostgreSQL |

SQLite works well for most individual users and does not require installing any additional software. Switch to PostgreSQL when you need concurrent write access from multiple processes or users, or when running a persistent server that handles multiple simultaneous requests.

## Verifying your configuration

Run `refgenie config get` to display the current configuration, including the active database backend:

```bash
refgenie config get
```

This shows the database URL (with credentials masked for PostgreSQL) along with other environment-based settings like `genome_folder` and `genome_stage_folder`.

## Related environment variables

| Variable | Default | Description |
|---|---|---|
| `REFGENIE_HOME_PATH` | `~/.refgenie` | Base directory for refgenie files |
| `REFGENIE_DB_CONFIG_PATH` | `$REFGENIE_HOME_PATH/refgenie_db_config.yaml` | Path to database config |
| `REFGENIE_GENOME_FOLDER` | `$REFGENIE_HOME_PATH/genomes` | Directory for genome assets |
| `REFGENIE_GENOME_STAGE_FOLDER` | `$REFGENIE_HOME_PATH/archives` | Directory for staged assets |
| `REFGENIE_LOG_LEVEL` | `INFO` | Log verbosity |
| `REFGENIE_REFGET_STORE_URL` | `None` | URL for an external RefgetStore |
