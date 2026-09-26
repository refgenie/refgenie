# Run a refgenie server

Refgenie includes a built-in server that exposes a REST API for listing, browsing, and downloading genome assets. The server also provides GA4GH DRS endpoints and can aggregate data channels from multiple sources.

## Start the server

```bash
refgenie serve
```

This starts the server on port 8000 by default.

### Flags

| Flag | Description |
|---|---|
| `-p`, `--port` | Port to run the server on (default: 8000) |
| `-r`, `--reload` | Enable auto-reload on code changes (for development) |

Example with custom port:

```bash
refgenie serve -p 5000
```

## Configuration via environment variables

The server reads the same environment variables as the CLI:

| Variable | Default | Description |
|---|---|---|
| `REFGENIE_HOME_PATH` | `~/.refgenie` | Base directory for refgenie files |
| `REFGENIE_DB_CONFIG_PATH` | `$REFGENIE_HOME_PATH/refgenie_db_config.yaml` | Path to database config |
| `REFGENIE_GENOME_FOLDER` | `$REFGENIE_HOME_PATH/genomes` | Directory for genome asset data |
| `REFGENIE_GENOME_STAGE_FOLDER` | `$REFGENIE_HOME_PATH/archives` | Directory for staged assets (served to clients) |
| `REFGENIE_LOG_LEVEL` | `INFO` | Log verbosity |

For production deployments, configure PostgreSQL as the database backend. See [Database backends](../database.md) for setup details.

## Staging workflow for serving

The server serves assets from the stage folder (`REFGENIE_GENOME_STAGE_FOLDER`). Before assets are available to clients, they must be built and staged:

```bash
# 1. Build the asset
refgenie build hg38/fasta

# 2. Stage the asset (creates symlinks or tarballs in the stage folder)
refgenie stage stage hg38/fasta

# 3. Start the server
refgenie serve
```

You can combine building and staging in one step:

```bash
refgenie build hg38/fasta --stage
```

See [Stage assets for serving](../staging.md) for details on the staging system.

## Interacting with the API

Once the server is running, you can:

- Browse the API docs at `http://localhost:8000/docs` (Swagger UI)
- List available assets via the REST API
- Use `refgenie pull` from any client machine pointed at your server

### Subscribe a client to a server

On a client machine, subscribe to a remote server:

```bash
refgenie subscribe -s http://your-server:8000
```

Then pull assets:

```bash
refgenie pull hg38/fasta
```

## Server vs. dashboard

| Feature | `refgenie serve` | `refgenie dash` |
|---|---|---|
| Purpose | Serve assets to remote clients via API | Browse your local refgenie instance |
| Network | Exposes HTTP API for external clients | Local web UI only |
| Default port | 8000 | 8080 |

See [Use the dashboard](../dash.md) for dashboard documentation.

## Related pages

- [Serving modes](../serving_modes.md) -- How asset classes declare their delivery strategy
- [Configure remote storage](../remotes.md) -- Push staged assets to S3 or other cloud storage
- [Database backends](../database.md) -- SQLite vs. PostgreSQL configuration
- [Use the dashboard](../dash.md) -- Local web UI for browsing genomes and assets
