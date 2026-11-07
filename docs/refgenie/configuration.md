# Configure refgenie

Refgenie keeps its settings in a database (SQLite or PostgreSQL) and reads a few environment variables for paths and behavior. The settings include the genome folder, the stage folder, subscribed servers, data channels, and remotes.

## Initialize

```bash
refgenie init
```

This creates `~/.refgenie/` with a SQLite database, a `genomes/` folder, and a database config file. It also subscribes to the default server, `https://api.refgenie.org`, so `refgenie pull` works right away. To start with no subscriptions, run `refgenie init --no-default-server` or set `REFGENIE_SUBSCRIBE_DEFAULT_SERVER=false`. Running `refgenie init` again leaves existing subscriptions as they are.

## Environment variables

| Variable | Default | Purpose |
|---|---|---|
| `REFGENIE_HOME_PATH` | `~/.refgenie` | Base directory for all refgenie files |
| `REFGENIE_LOG_LEVEL` | `INFO` | Logging verbosity |
| `REFGENIE_GENOME_FOLDER` | `$REFGENIE_HOME_PATH/genomes` | Directory for genome assets |
| `REFGENIE_GENOME_STAGE_FOLDER` | `$REFGENIE_HOME_PATH/archives` | Directory for staged assets |
| `REFGENIE_DB_CONFIG_PATH` | `$REFGENIE_HOME_PATH/refgenie_db_config.yaml` | Path to the database config file |
| `REFGENIE_SUBSCRIBE_DEFAULT_SERVER` | `true` | Subscribe a new configuration to `https://api.refgenie.org` |
| `REFGENIE_CONTAINER_ENGINE` | `auto` | What runs a recipe's container image: `docker`, `apptainer`, or `auto` (Docker when its daemon answers, else Apptainer) |
| `REFGENIE_APPTAINER_CACHEDIR` | unset | Passed to Apptainer as `APPTAINER_CACHEDIR`, so image pulls land there instead of under `$HOME` |

For example:

```bash
export REFGENIE_GENOME_FOLDER=~/my_genomes
export REFGENIE_DB_CONFIG_PATH=~/my_refgenie_db_config.yaml
```

Coming from refgenie 0.x? The `$REFGENIE` variable that pointed to the 0.x config file is ignored in 1.x, which uses `REFGENIE_HOME_PATH` and `REFGENIE_DB_CONFIG_PATH` instead. To move an existing setup, see [Upgrade from legacy](upgrade.md).

The [README](https://github.com/refgenie/refgenie#environment-variables) lists the variables for the server, dashboard, and plugins too.

## Database

The database config file selects SQLite or PostgreSQL. See [Database backends](database.md) for the file format and when to use each.

## Servers

Subscribe to servers to pull assets and look up remote paths. You can also manage subscriptions on the **Servers** page of `refgenie dash`.

```bash
refgenie subscribe -s https://my-lab-server.example.org
refgenie unsubscribe -s https://my-lab-server.example.org
```

## Inspect the configuration

```bash
refgenie config get
```

This prints a table of the stored settings (genome folder, stage folder, subscribed servers, plugin settings), then the settings that come from environment variables:

```
log_level: LogLevel.INFO
genome_folder: /home/user/.refgenie/genomes
genome_stage_folder: /home/user/.refgenie/archives
database_config_path: /home/user/.refgenie/refgenie_db_config.yaml
...
```

## Configure from Python

The `Refgenie` constructor takes these options:

- `database_config_path`: path to a database config file, instead of `REFGENIE_DB_CONFIG_PATH`.
- `database_engine`: a SQLAlchemy engine you already have, such as a pooled connection or a test database. If you pass both, the engine wins.
- `server_clients_mapping`: a dict of server URL to client object, for custom or mock server clients. Each client is a `RefgenieserverClient` from `refgenie.managers.sources.client`.

```python
from sqlalchemy import create_engine
from refgenie import Refgenie

engine = create_engine("postgresql+psycopg://user:password@host:5432/dbname")
rgc = Refgenie(database_engine=engine)
```
