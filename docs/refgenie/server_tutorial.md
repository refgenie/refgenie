# Running a Refgenie Server

This tutorial walks you through setting up and running your own refgenie server. A refgenie server allows you to share genome assets with others on your network or the public internet.

## Why Run Your Own Server?

There are several reasons you might want to run your own refgenie server:

- **Internal distribution**: Share assets privately across your organization without uploading to public servers
- **Custom assets**: Serve specialized or proprietary genome assets that aren't available on public servers
- **Performance**: Run a local server for faster asset access within your network
- **Community sharing**: Distribute your resources through the familiar refgenie interface

## Prerequisites

Before starting a server, you need:

1. **Refgenie installed with server extras**:
   ```bash
   pip install refgenie[server]
   ```

2. **Built assets**: You need local assets to serve. If you haven't built any yet, see the [build documentation](build.md).

3. **Staged assets**: Assets must be staged before they can be served. See [Stage assets for serving](staging.md).

!!! info "Server vs Dashboard"
    Refgenie provides two web interfaces:

    - **`refgenie serve`** - A full server for distributing assets to remote clients
    - **`refgenie dash`** - A local web UI for browsing and managing your own assets (see [Use the dashboard](dash.md))

## Staging assets

Before serving assets, you must stage them. Staging puts a servable copy of each asset (a symlink, a `.tgz` archive, or both, depending on the asset class) into the stage folder.

### Set up a stage folder

The stage folder defaults to `~/.refgenie/archives`. You can change it during initialization or with an environment variable:

```bash
# Set via environment variable
export REFGENIE_GENOME_STAGE_FOLDER=~/.refgenie/archives

# Or specify during init
refgenie init --genome-stage-folder ~/.refgenie/archives
```

### Stage assets

Stage one or more built assets:

```bash
refgenie stage add hg38/fasta hg38/bowtie2_index mm10/fasta
```

### Stage during build

You can also stage an asset automatically when you build it:

```bash
refgenie build hg38/fasta --stage
```

### List staged assets

```bash
refgenie stage list
```

!!! warning "Staging Required for Serving"
    The server only serves staged assets. If an asset isn't staged, clients won't be able to pull it.

## Starting the Server

### Basic server startup

Start the server with default settings (port 8000):

```bash
refgenie serve
```

The server will be available at `http://localhost:8000`.

### Custom port

Specify a different port:

```bash
refgenie serve --port 8080
```

### Development mode with auto-reload

For development, enable auto-reload to automatically restart when code changes:

```bash
refgenie serve --reload
```

!!! tip "Production Deployment"
    For production use, consider running behind a reverse proxy (nginx, Apache) with HTTPS enabled. The `--reload` flag should not be used in production.

## Using the Local Dashboard

The dashboard is a local web interface for browsing and managing your own assets:

```bash
refgenie dash
```

By default, this opens `http://localhost:8080` in your browser. Specify a different port with:

```bash
refgenie dash --port 9000
```

See [Use the dashboard](dash.md) for what it can do.

## Server API Overview

The refgenie server exposes a REST API with several endpoint groups:

### Core Endpoints

| Endpoint | Description |
|----------|-------------|
| `/` or `/index` | Landing page with all available genomes |
| `/healthcheck` | Health check endpoint (returns `{"status": "ok"}`) |
| `/summary` | Database summary with counts of genomes, asset groups, and assets |
| `/species/summary` | Statistics grouped by species |

### Asset Endpoints

| Endpoint | Description |
|----------|-------------|
| `/page/genome/{digest}` | Genome splash page with available assets |
| `/page/asset/{digest}/{asset_group}` | Asset detail page |
| `/archives/{digest}/download` | Download an archive by its digest |

### GA4GH DRS Endpoints

The server implements the [GA4GH Data Repository Service (DRS)](https://ga4gh.github.io/data-repository-service-schemas/) specification:

| Endpoint | Description |
|----------|-------------|
| `/ga4gh/drs/v1/service-info` | DRS service information |
| `/ga4gh/drs/v1/objects/{object_id}` | Get DRS object metadata |
| `/ga4gh/drs/v1/objects/{object_id}/access/{access_id}` | Get access URL for download |

### Data Channel Endpoints

| Endpoint | Description |
|----------|-------------|
| `/data_channel/` | List available data channels |
| `/data_channel/index.yaml` | Compiled index of all data channels |

### Interactive API Documentation

Visit `/docs` for Swagger UI or `/redoc` for ReDoc-style documentation.

## Configuration for Production

### Using PostgreSQL

For production deployments, consider using PostgreSQL instead of SQLite. Create a database configuration file:

```yaml
# ~/.refgenie/refgenie_db_config.yaml
type: postgresql
name: refgenie
host: localhost
port: 5432
user: refgenie_user
password: your_secure_password
```

Set the path to this config:

```bash
export REFGENIE_DB_CONFIG_PATH=~/.refgenie/refgenie_db_config.yaml
```

See the [configuration documentation](configuration.md) for more details.

### Environment Variables

Key environment variables for server configuration:

| Variable | Default | Description |
|----------|---------|-------------|
| `REFGENIE_HOME_PATH` | `~/.refgenie` | Base directory for refgenie files |
| `REFGENIE_GENOME_FOLDER` | `$REFGENIE_HOME_PATH/genomes` | Directory for genome assets |
| `REFGENIE_GENOME_STAGE_FOLDER` | `$REFGENIE_HOME_PATH/archives` | Directory for staged assets |
| `REFGENIE_DB_CONFIG_PATH` | `$REFGENIE_HOME_PATH/refgenie_db_config.yaml` | Path to database config |
| `REFGENIE_LOG_LEVEL` | `INFO` | Logging verbosity |

### Serving from Remote Storage

For serving assets from cloud storage (e.g., S3), configure a remote:

```bash
refgenie remote add --type http --prefix https://your-bucket.s3.amazonaws.com/archives --name s3-archive
```

When a remote is configured, the server will redirect download requests to the remote URL instead of serving files directly.

### Configuring Data Channels

Data channels allow you to aggregate asset classes and recipes from external sources. Set the path to your data channels configuration:

```bash
export DATA_CHANNELS_CONFIG_PATH=/path/to/data_channels.yaml
```

See [data channels documentation](data_channels.md) for the configuration format.

## Docker Deployment

For containerized deployments, you can run the server with Docker:

```bash
docker run --rm -d -p 80:80 \
    -v /path/to/genomes:/genomes \
    -v /path/to/archives:/archives \
    -e REFGENIE_GENOME_FOLDER=/genomes \
    -e REFGENIE_GENOME_STAGE_FOLDER=/archives \
    --name refgenie-server \
    refgenie/refgenie:latest refgenie serve --port 80
```

## Client Configuration

Once your server is running, clients can subscribe to it:

```bash
refgenie subscribe http://your-server-address:8000
```

Then pull assets normally:

```bash
refgenie pull hg38/fasta
```

## Troubleshooting

### Server won't start

- Ensure the `server` extras are installed: `pip install refgenie[server]`
- Check that the port isn't already in use
- Verify database configuration is correct

### Assets not showing up

- Confirm assets are staged: `refgenie stage list`
- Check the stage folder path is correct
- Verify database contains the assets: `refgenie list`

### Download failures

- Ensure staged files exist in the stage folder
- If using remote storage, verify the remote configuration is correct
- Check server logs for detailed error messages

## Next Steps

- [Configuration options](configuration.md) - Detailed configuration reference
- [Data channels](data_channels.md) - Set up external data sources
- [Building assets](build.md) - Create assets to serve
- [Public servers](servers.md) - List of available public servers
