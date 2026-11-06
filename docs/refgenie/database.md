# Database backends

Refgenie keeps asset metadata, genome records and settings in a database: SQLite (the default) or PostgreSQL.

## SQLite (default)

SQLite needs no setup. The first time refgenie runs, it writes a config file that points at a SQLite database in `~/.refgenie/refgenie`:

```yaml
# ~/.refgenie/refgenie_db_config.yaml
type: sqlite
path: ~/.refgenie/refgenie
```

## PostgreSQL

Use PostgreSQL for shared servers and for anything that writes from several processes at once. Write the config file yourself:

```yaml
# ~/.refgenie/refgenie_db_config.yaml
type: postgresql
name: refgenie
user: refgenie_user
password: your_secure_password
host: localhost
port: 5432
```

Values can use environment variables, such as `password: $POSTGRES_PASSWORD`. Refgenie connects with the `psycopg` driver.

## Choosing a backend

| Scenario | Backend |
|---|---|
| One user on a workstation | SQLite |
| CI or a short-lived environment | SQLite |
| A shared lab or team database | PostgreSQL |
| A refgenie server (`refgenie serve`) | PostgreSQL |
| A cloud deployment with several workers | PostgreSQL |

## Config file location

Refgenie reads the config file from `REFGENIE_DB_CONFIG_PATH`, which defaults to `$REFGENIE_HOME_PATH/refgenie_db_config.yaml`:

```bash
export REFGENIE_DB_CONFIG_PATH=/path/to/refgenie_db_config.yaml
```

`refgenie config get` shows the settings refgenie is using, including this path, `genome_folder` and `genome_stage_folder`.

## Related environment variables

| Variable | Default | Description |
|---|---|---|
| `REFGENIE_HOME_PATH` | `~/.refgenie` | Base directory for refgenie files |
| `REFGENIE_DB_CONFIG_PATH` | `$REFGENIE_HOME_PATH/refgenie_db_config.yaml` | Path to the database config |
| `REFGENIE_GENOME_FOLDER` | `$REFGENIE_HOME_PATH/genomes` | Directory for genome assets |
| `REFGENIE_GENOME_STAGE_FOLDER` | `$REFGENIE_HOME_PATH/archives` | Directory for staged assets |
| `REFGENIE_LOG_LEVEL` | `INFO` | Log verbosity |
| `REFGENIE_STORES` | unset | Server only: a JSON list of RefgetStores to register at startup. See [`deployment/stores.example.json`](https://github.com/refgenie/refgenie/blob/master/deployment/stores.example.json). |
