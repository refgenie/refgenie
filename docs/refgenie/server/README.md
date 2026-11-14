# Server reference

`refgenie serve` runs the refgenie server: a REST API for listing, browsing and downloading assets, plus GA4GH DRS and sequence collections endpoints. For a step-by-step setup, see [Run a server](../server_tutorial.md).

## Command

```text
refgenie serve [-p PORT] [-r]
```

| Flag | Description |
|---|---|
| `-p`, `--port` | Port to listen on (default `8000`). The server listens on all interfaces. |
| `-r`, `--reload` | Restart when the code changes. For development only. |

It needs the `server` extra: `pip install "refgenie[server]>=1.0.0a2"`.

## Endpoints

The full list, with request and response formats, is at `/docs` (Swagger UI) and `/redoc` on any running server.

| Path | What it serves |
|---|---|
| `/` | A short landing page with links to the docs |
| `/v4/genomes`, `/v4/assets`, `/v4/asset_groups`, `/v4/asset_classes`, `/v4/recipes`, `/v4/aliases` | Catalog listings, with [pagination and search](pagination_and_search.md) |
| `/v4/archives/{asset_digest}/download` | An archive-mode asset as a `.tgz` |
| `/v4/assets/{asset_digest}/files/{file_path}` | One file of a file-mode asset |
| `/v4/summary`, `/v4/species/summary` | Counts of genomes, asset groups and assets |
| `/v4/data_channel/` | The data channels this server re-serves (see [Serve data channels](data_channels.md)) |
| `/ga4gh/drs/v1/` | [GA4GH Data Repository Service](https://ga4gh.github.io/data-repository-service-schemas/) objects for staged assets |
| `/seqcol` | [GA4GH sequence collections](https://ga4gh.github.io/refget/seqcols/) for the server's genomes |
| `/service-info` | What this server is and which features it has |
| `/mcp` | A read-only [MCP server](../mcp.md) for AI assistants |
| `/v4/healthcheck` | Returns `{"status": "ok"}` |

Downloads redirect to a remote when the asset was pushed to one. See [Remotes on a server](remotes.md).

## Environment variables

The server reads the same variables as the CLI, plus a few of its own.

| Variable | Default | Description |
|---|---|---|
| `REFGENIE_HOME_PATH` | `~/.refgenie` | Base directory for refgenie files |
| `REFGENIE_DB_CONFIG_PATH` | `$REFGENIE_HOME_PATH/refgenie_db_config.yaml` | Database config. Use PostgreSQL for a production server; see [Database backends](../database.md). |
| `REFGENIE_GENOME_FOLDER` | `$REFGENIE_HOME_PATH/genomes` | Where asset data lives |
| `REFGENIE_GENOME_STAGE_FOLDER` | `$REFGENIE_HOME_PATH/archives` | Where staged assets live |
| `REFGENIE_LOG_LEVEL` | `WARNING` | Log verbosity |
| `REFGENIE_COLLECTIONS` | unset | A JSON list of collection pages for the web UI, each drawn from one or more declared stores. Entries naming an unknown store are skipped. See [`deployment/collections.example.json`](https://github.com/refgenie/refgenie/blob/master/deployment/collections.example.json). |
| `REFGENIE_CATALOG_URL` | unset | A publish catalog (made by `refgenie catalog-export`) to import at startup and then on a schedule. The catalog also lists the RefgetStores the server reads; a change to that list takes effect at the next import, without a restart. |
| `CATALOG_IMPORT_INTERVAL_SECONDS` | `86400` | How often to re-import `REFGENIE_CATALOG_URL` |
| `DATA_CHANNELS_CONFIG_PATH` | the packaged `data_channels.yaml` | Which data channels to re-serve |
| `DOWNLOAD_COUNT_DUMP_JOB_INTERVAL_SECONDS` | `60` | How often download counts are saved |
| `REFGENIE_SERVER_PLUGINS` | `false` | Run [plugins](../plugins.md) on the server |

## Download counts

The server counts archive downloads (`/v4/archives/{asset_digest}/download`) in memory and adds them to the database every `DOWNLOAD_COUNT_DUMP_JOB_INTERVAL_SECONDS`. The totals appear as `download_count` in `/v4/archives`.

## Species names (`species.json`)

A store may carry an optional `species.json` next to its `rgstore.json`. The server uses it to show common names and to let users search genomes by any species name.

```json
{"species": [
  {"taxon_id": 10116, "scientific_name": "Rattus norvegicus", "common_name": "Norway rat", "synonyms": ["brown rat", "rat"]}
]}
```

- `taxon_id`: NCBI taxon id, an integer, required, unique in the file.
- `scientific_name`: required.
- `common_name` and `synonyms` (a list of other names): optional.

`refgenie store sync` reads the file from every enabled store, local or http(s). A missing file is fine and bad entries are skipped. If a store's file cannot be read, the existing species rows are kept.

Stores are merged in priority order. The scientific name comes from the highest-priority store, the common name is the first non-blank one, and every other name from every store is kept once (ignoring case). Each listed taxon is rewritten, so a name removed from a file disappears. Taxa that no file lists are kept; delete them by hand if needed. A server's catalog import updates species the same way.

## Related pages

- [Run a server](../server_tutorial.md)
- [Serving modes](../serving_modes.md): how asset classes choose file or archive delivery
- [Stage assets for serving](../staging.md)
- [Use the dashboard](../dash.md): the local web UI for your own refgenie
