# Refgenie Configuration

Refgenie uses a flexible configuration system to manage reference genome assets, supporting both local and remote operation modes. This document describes how to configure refgenie, including environment variables, database backends, and advanced options.

## Configuration Overview

Refgenie stores its configuration in a database (SQLite or PostgreSQL) and uses a set of environment variables to control paths and behavior. The configuration tracks:

- The genome folder (where assets are stored)
- The genome archive folder (for asset archives)
- Database connection settings
- Subscribed remote servers
- Registered data channels and remotes

A visual overview of the configuration and database model is available in the database model diagram.

## Initialization

To initialize refgenie with default settings, simply run:

```bash
refgenie1 init
```

This will create a default configuration in your home directory under `.refgenie/` and set up a local SQLite database. You can override the default locations and database settings using environment variables (see below).

## Environment Variables

Refgenie supports the following environment variables for configuration:

| Variable                      | Default                                 | Purpose                                                      |
|-------------------------------|-----------------------------------------|--------------------------------------------------------------|
| `REFGENIE_HOME_PATH`          | `~/.refgenie`                           | Base directory for all refgenie files                        |
| `REFGENIE_LOG_LEVEL`          | `INFO`                                  | Logging verbosity                                            |
| `REFGENIE_GENOME_FOLDER`      | `$REFGENIE_HOME_PATH/genomes`           | Directory for genome assets                                  |
| `REFGENIE_GENOME_ARCHIVE_FOLDER` | `$REFGENIE_HOME_PATH/archives`        | Directory for asset archives                                 |
| `REFGENIE_DB_CONFIG_PATH`     | `$REFGENIE_HOME_PATH/refgenie_db_config.yaml` | Path to the database config YAML file (see next section to learn more about the schema) |

You can set these variables in your shell before running refgenie, for example:

```zsh
export REFGENIE_GENOME_FOLDER=~/my_genomes
export REFGENIE_DB_CONFIG_PATH=~/my_refgenie_db_config.yaml
```

## Database Backend

Refgenie supports both SQLite (default) and PostgreSQL for storing metadata. The database connection is configured via a YAML file, whose path is set by `REFGENIE_DB_CONFIG_PATH`.

### Database Configuration YAML Schema

The database backend is configured via a YAML file whose path is set by `REFGENIE_DB_CONFIG_PATH`. The schema and required fields depend on the database type you choose. Refgenie uses a **discriminated union** to select the correct configuration model based on the `type` field in the YAML file:

- If `type: sqlite`, the SQLite configuration schema is used.
- If `type: postgresql`, the PostgreSQL configuration schema is used.

This means the `type` field acts as a discriminator, and the rest of the fields are validated according to the selected backend.

#### Common Fields

- `type`: The database backend type. Must be either `sqlite` or `postgresql`. This is the discriminator for the schema.

#### SQLite Example

```yaml
type: sqlite
path: ~/refgenie_db/refgenie  # Path to the SQLite database file
```

- Only `type` and `path` are required for SQLite.
- `path` should be a valid file path where the SQLite database will be stored.

#### PostgreSQL Example

```yaml
type: postgresql
name: refgenie                # Name of the PostgreSQL database
host: localhost               # Hostname of the PostgreSQL server
port: 5432                    # Port number
user: postgres                # Username
password: mysecretpassword    # Password
```

- All fields shown above are required for PostgreSQL.
- `name` is the database name (not a file path).
- `host`, `port`, `user`, and `password` must be specified to connect to the server.

> **Note:** Refgenie will create a default config file if one does not exist. You can edit this file to point to your preferred backend and provide the necessary connection details. The config loader will automatically select the correct schema based on the `type` field using a discriminated union.


### Why use a non-local PostgreSQL database?

While SQLite is simple and works well for single-user or small-scale setups, configuring refgenie to use a PostgreSQL database (especially a remote or managed instance) provides several advantages for larger or more demanding environments:

- **High performance for concurrent access**: PostgreSQL supports multiple simultaneous users and high-throughput operations, making it suitable for large teams or automated pipelines.
- **Distributed and multi-user workflows**: A remote PostgreSQL server allows multiple users, compute nodes, or cloud services to access and update the same refgenie database, supporting collaborative and distributed workflows.
- **Scalability for large organizations**: PostgreSQL can efficiently handle large datasets and many assets, making it ideal for institutional or enterprise deployments.
- **High availability and reliability**: Managed PostgreSQL services (e.g., AWS RDS, Google Cloud SQL) offer automated backups, failover, and monitoring, reducing downtime and risk of data loss.
- **Centralized management**: A single remote database can serve as the authoritative source of truth for all reference assets across an organization.

For most individual users or small labs, SQLite is sufficient. For production, cloud, or institutional deployments, PostgreSQL is recommended.

## Local and Remote Configuration

- **Local mode**: All assets and metadata are managed on your local machine using the configured database and folders.
- **Remote mode**: You can subscribe to remote refgenie servers to pull assets or use remote seek functionality. Use the CLI or Python API to add server subscriptions:

```bash
refgenie1 subscribe http://refgenomes.databio.org
refgenie1 unsubscribe -s http://refgenomes.databio.org
```

You can list, add, or remove remote servers and data channels at any time.

## Inspecting Configuration

You can inspect the current configuration using:

```bash
refgenie1 config get
```

This will print a table with the current settings, including environment-based overrides.

```
                                                   Refgenie configuration
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ genome_folder                       ┃ version ┃ genome_archive_folder                ┃ servers                           ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ /Users/johndoe/.refgenie/genomes    │ 1       │ /Users/johndoe/.refgenie/archives    │ ['http://refgenomes.databio.org'] │
└─────────────────────────────────────┴─────────┴──────────────────────────────────────┴───────────────────────────────────┘

Environment-based configuration

log_level: LogLevel.INFO
genome_folder: /Users/johndoe/.refgenie/genomes
genome_archive_folder: /Users/johndoe/.refgenie/archives
database_config_path: /Users/johndoe/.refgenie/refgenie_db_config.yaml
```

## Advanced configuration

### Data Channels and Remotes

Refgenie supports registration of data channels (external sources of asset classes and recipes) and remotes (external storage backends, e.g., S3). These are tracked in the configuration database and can be managed via CLI or API.

### Advanced: Programmatic Configuration and Customization

When using refgenie as a Python library, users can inject custom objects into the `Refgenie` class initializer to control its behavior beyond environment variables and config files.

### Provide path to database config file

When using refgenie as a Python library, you can specify the path to the database configuration file directly in the `Refgenie` constructor with `database_config_path`. This allows you to programmatically set up the database connection without relying on environment variables.

#### Custom Database Engine

You can provide a custom SQLAlchemy database engine to the `Refgenie` constructor via the `database_engine` argument. This allows you to:

- Use a pre-configured or pooled database connection
- Integrate with existing infrastructure or testing environments

Example:

```python
from refgenie import Refgenie
from sqlalchemy import create_engine

engine = create_engine("postgresql://user:password@host:port/dbname")
refgenie = Refgenie(database_engine=engine)
```

> Note that this argument and the `database_config_path` argument are mutually exclusive. If both are provided, `database_config_path` will be ignored.

#### Custom Server Clients

You can also inject a mapping of server URLs to custom server client objects using the `server_clients_mapping` argument. Importantly, each client must implement the `ServerClient` protocol (see `refgenie.server.models.ServerClient`). 

This is useful for:

- Extending or **customizing remote server interactions**
- Using mock or test clients

Example:

```python
from refgenie import Refgenie
from mymodule import MyCustomServerClient

server_clients = {
    "http://myserver": MyCustomServerClient(),
}
refgenie = Refgenie(server_clients_mapping=server_clients)
```
```
╭──────────────────────────────────────── ServerClient Protocol structure ────────────────────────────────────────╮
│ class ServerClient(*args, **kwargs):                                                                            │
│                                                                                                                 │
│ Protocol for the server client.                                                                                 │
│                                                                                                                 │
│             server_url = <property object at 0x11b0bf880>                                                       │
│ download_with_progress = def download_with_progress(self, operation_id: str, output_path: pathlib.Path, params: │
│                          Optional[Dict] = None, url_format_params: Optional[Dict[str, str]] = None, name:       │
│                          Optional[str] = None) -> pathlib.Path: Download asset at given URL to given filepath,  │
│                          show progress along the way.                                                           │
│                    get = def get(self, operation_id: str, params: Optional[Dict] = None, url_format_params:     │
│                          Optional[Dict[str, str]] = None) -> Union[Dict, str]: Send a GET request to the        │
│                          specified operation ID.                                                                │
╰─────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯
```

> **Note:** If a provided client does not implement the required protocol, refgenie will raise a `ValueError`.

This programmatic approach enables full control for advanced, production, or testing scenarios.

### Database Migrations

Refgenie uses Alembic for database schema migrations. For advanced users or developers, see the migrations documentation in the [refgenie repository](https://github.com/refgenie/refgenie1) for details on managing schema changes.
