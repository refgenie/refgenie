# Refgenie 1.0

[![Run tests](https://github.com/refgenie/refgenie/actions/workflows/test.yaml/badge.svg)](https://github.com/refgenie/refgenie/actions/workflows/test.yaml)
[![PyPI](https://img.shields.io/pypi/v/refgenie.svg)](https://pypi.org/project/refgenie/)
[![Python versions](https://img.shields.io/pypi/pyversions/refgenie.svg)](https://pypi.org/project/refgenie/)
[![License: BSD-2-Clause](https://img.shields.io/badge/license-BSD--2--Clause-blue.svg)](https://github.com/refgenie/refgenie/blob/master/LICENSE.txt)

> ## ⚠️ This is an alpha release
>
> Refgenie 1.0 is **pre-release software**. The CLI, the Python API, the REST
> API, and the database schema can all change without a deprecation period, and
> a future alpha may require you to rebuild your assets. Do not depend on it in
> production pipelines yet. Please do try it and
> [report what breaks](https://github.com/refgenie/refgenie/issues).

Refgenie is a reference genome asset manager. It organizes, builds, and serves
the reference genome files that bioinformatics workflows depend on — FASTA
files, aligner indexes, annotations — and gives every one of them a stable
identifier you can resolve from a script.

## Features

- **Modular asset system**: manage different types of reference genome assets (FASTA, BWA index, etc.)
- **Content-addressed assets**: an asset's digest is derived from the files it contains
- **Seqcol genome identity**: genomes are identified by GA4GH refget sequence-collection digests
- **Database backend**: SQLite or PostgreSQL for asset metadata, with alembic migrations
- **Data channels**: register and sync asset classes and recipes from external sources
- **Server and web UI**: a REST API (with GA4GH DRS endpoints) and a local web UI
- **MCP server**: expose refgenie to LLM agents over the Model Context Protocol
- **Self-describing to agents**: every instance serves `/SKILL.md`, a capability
  doc that teaches an assistant this instance's HTTP API
- **Command-line interface**: build, seek, pull, push, and manage genome assets

## Installation

Install from [PyPI](https://pypi.org/project/refgenie/):

```bash
pip install refgenie
```

Because this is an alpha, `pip` needs to be told to accept a pre-release
version:

```bash
pip install --pre refgenie
```

### Optional extras

The base install gives you the CLI and the Python API. Everything else is an
extra:

| Extra | Install | What it enables | Adds |
| --- | --- | --- | --- |
| *(none)* | `pip install refgenie` | CLI, Python API, building, pulling, staging | — |
| `dash` | `pip install "refgenie[dash]"` | `refgenie dash`, the local web UI | fastapi, uvicorn |
| `server` | `pip install "refgenie[server]"` | `refgenie serve`, the public REST API + DRS + `/mcp` endpoint (includes `dash`) | apscheduler, mcp |
| `mcp` | `pip install "refgenie[mcp]"` | `refgenie-mcp`, the stdio MCP server for LLM agents | mcp |
| `snakemake` | `pip install "refgenie[snakemake]"` | snakemake-based bulk asset building | snakemake |

`refgenie serve`, `refgenie dash`, and `refgenie-mcp` all appear in
`refgenie --help` on a base install. Running one without its extra prints a
message telling you which extra to install.

PostgreSQL support needs no extra — `psycopg` is a base dependency.

## Quick Start

Run these in order — refgenie ships with no recipes, so the data channel steps
are required before anything can be built.

```bash
# 1. Create the refgenie home, database, and genome folders
refgenie init

# 2. Register and sync a data channel; this is where asset classes and
#    recipes (including the 'fasta' recipe) come from
refgenie data-channel add refgenie https https://refgenie.github.io/refgenie-registry/index.yaml
refgenie data-channel sync refgenie --exists-ok

# 3. Register a genome from a FASTA file. Because the 'fasta' recipe is now
#    available, this also builds the genome's fasta asset.
refgenie genome init --fasta genome.fa.gz --name genome_name

# 4. Build any other asset for that genome (idempotent; skips if it exists)
refgenie build genome_name/fasta --files fasta=genome.fa.gz

# 5. Inspect what you have, and get a path to it
refgenie list
refgenie seek genome_name/fasta
```

`refgenie list` prints something like:

```
                      Refgenie assets. Source: local
┏━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━┓
┃ Aliases     ┃ Genome digest                    ┃ Asset group ┃ Asset   ┃
┡━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━┩
│ genome_name │ jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP │ fasta       │ default │
└─────────────┴──────────────────────────────────┴─────────────┴─────────┘
```

## Environment variables

Refgenie1 reads the following environment variables. All are optional — defaults
under `~/.refgenie/` are used when unset.

| Variable | Default | Purpose |
| --- | --- | --- |
| `REFGENIE_HOME_PATH` | `~/.refgenie` | Root directory for refgenie1 state. Auto-created at first import. |
| `REFGENIE_DB_CONFIG_PATH` | `$REFGENIE_HOME_PATH/refgenie_db_config.yaml` | Path to the refgenie1 database connection YAML. |
| `REFGENIE_GENOME_FOLDER` | `$REFGENIE_HOME_PATH/genomes` | Where built asset files live. Auto-created on `refgenie init`. |
| `REFGENIE_GENOME_STAGE_FOLDER` | `$REFGENIE_HOME_PATH/archives` | Where staged/archived assets live before push. Auto-created on `refgenie init`. |
| `REFGENIE_STORES` | unset | JSON array of federated refget stores registered on server boot, e.g. `[{"name": "jungle", "url": "https://.../jungle/", "priority": 10}]`. Lower priority wins alias ties. See `deployment/stores.example.json`. |
| `REFGENIE_BRIDGE_MODE` | `read` | Localhost bridge policy for `refgenie dash`. `off` = no cross-origin access; `read` = allowlisted public origins may read `/ping`, `/v4`, `/v1/jobs`; `full` = additionally allows `POST /v1/actions/pull`. `refgenie dash --bridge` overrides it for one run. |
| `REFGENIE_BRIDGE_ORIGINS` | `https://refgenie.org` | Comma-separated exact origins the bridge admits. **Replaces** the default rather than appending. |
| `REFGENIE_BRIDGE_ORIGIN_REGEX` | unset | Origin regex, for SPA dev servers and preview deployments only. Wildcard-equivalent patterns are refused at startup. |
| `REFGENIE_BRIDGE_EXPOSE_PATHS` | `false` | Whether `/ping` may report the genome folder path. Off by default: a home-directory path leaks the OS username to every allowlisted origin. |
| `REFGENIE_WEB_DIST` | unset | Override the built web UI directory. Default is the packaged `refgenie/server/webui/`. |
| `REFGENIE_LOG_LEVEL` | `INFO` | Python logging level. |
| `REFGENIE_ENCRYPTION_KEY` | (auto-generated, dev only) | Fernet key for encrypting credentials. **Set explicitly in production.** |
| `REFGENIE_DISABLE_PLUGINS` | unset | `1`/`true`/`all` turns off every plugin; a comma-separated list of plugin names turns off those. Wins over everything else. |
| `REFGENIE_SERVER_PLUGINS` | `false` | Run plugins in server mode (`refgenie serve`). Off by default; local mode and `refgenie dash` always run them. |

Migrating from refgenie 0.x? Most users should run the `refgenie-upgrade` tool -- see the [upgrade walkthrough](https://docs.refgenie.org/refgenie/upgrade/). For a reference of all API and CLI breaking changes, see the [migration guide](https://github.com/refgenie/refgenie/blob/master/docs/refgenie/migration-from-legacy.md).

## Usage

For a complete list of commands:

```bash
refgenie --help
```

Longer walkthroughs:

- [CLI documentation](https://github.com/refgenie/refgenie/blob/master/docs/refgenie/refgenie-cli.md) (work in progress)
- [Tutorial notebook](https://docs.refgenie.org/refgenie/notebooks/refgenie/)

## Data Channels

Refgenie ships with no recipes or asset classes. A **data channel** is an index
of them that you register and sync; the canonical one is published by
[refgenie-registry](https://github.com/refgenie/refgenie-registry).

```bash
# Add a data channel: <name> <type> <index url>
refgenie data-channel add refgenie https https://refgenie.github.io/refgenie-registry/index.yaml

# Sync its asset classes and recipes into your local database
refgenie data-channel sync refgenie --exists-ok

# See what you have
refgenie recipe list
refgenie asset-class list
```

> Recipes contain commands, and building an asset runs them. Only sync channels
> you trust.

## Server, dashboard, and MCP

```bash
# Local web UI (needs the 'dash' extra)
refgenie dash

# Public REST API, including GA4GH DRS endpoints and an /mcp endpoint
# (needs the 'server' extra)
refgenie serve --port 8000

# stdio MCP server for LLM agents (needs the 'mcp' or 'server' extra)
refgenie-mcp
```

With the server running, interactive API docs are at `/docs` and the OpenAPI
document is at `/openapi.json`. `/SKILL.md` is the same API written for an AI
assistant; it ships in the web UI bundle (`frontend/public/SKILL.md`), so both
modes serve it. The HTTP MCP endpoint is at `/mcp/mcp` — the sub-application is
mounted at `/mcp` and registers its own `/mcp` route inside itself.

`refgenie dash` and `refgenie serve` are the same application in two modes: the
same JSON API under `/v4`, the same web UI, and a local-only command surface
under `/v1` that only the dash has. The UI asks `/service-info` which mode it is
talking to and shows only what that mode can do. `refgenie dash` binds
`127.0.0.1` only; use `refgenie serve` for anything reachable over a network.

A wheel ships the built web UI. In a source checkout, build it once with
`npm --prefix frontend run build` — without it the API still works and `/`
explains what to run.

## Development

This project uses [uv](https://github.com/astral-sh/uv) for package management and [Task](https://taskfile.dev/) for development tasks. See [docs/refgenie/development.md](https://github.com/refgenie/refgenie/blob/master/docs/refgenie/development.md) for the full developer quickstart.

```bash
# Clone and install
git clone https://github.com/refgenie/refgenie.git
cd refgenie1
uv sync

# Run tests
task test

# Format code
task reformat

# Run refgenie CLI
uv run refgenie <command>
```

## Contributing

Contributions are welcome. Please read
[CONTRIBUTING.md](https://github.com/refgenie/refgenie/blob/master/CONTRIBUTING.md)
for branch conventions, what a reviewable pull request looks like, and how to
run both test suites, and
[docs/refgenie/development.md](https://github.com/refgenie/refgenie/blob/master/docs/refgenie/development.md)
for the environment setup.

## License

Refgenie is released under the BSD 2-Clause License. See
[LICENSE.txt](https://github.com/refgenie/refgenie/blob/master/LICENSE.txt).

## Links

- [Source code](https://github.com/refgenie/refgenie)
- [Issue tracker](https://github.com/refgenie/refgenie/issues)
- [Changelog](https://github.com/refgenie/refgenie/blob/master/CHANGELOG.md)
- [Documentation](https://github.com/refgenie/refgenie/tree/master/docs)
- [Data channel registry](https://github.com/refgenie/refgenie-registry)
- [Migration guide](https://github.com/refgenie/refgenie/blob/master/docs/refgenie/migration-from-legacy.md) (for refgenie 0.x users)