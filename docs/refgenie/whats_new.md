# What's new in refgenie 1.0

Refgenie 1.0 is a ground-up rewrite that replaces the three legacy packages (refgenie, refgenconf, refgenieserver) with a single unified package. It installs as a new package rather than an in-place upgrade of the legacy tools, but existing legacy (v0.3/v0.4) installs can be migrated automatically -- see [Upgrading from legacy](upgrade.md).

## What's new for users

### Custom recipes and assets

The set of asset types and the recipes for building them were previously hardcoded in refgenie -- adding a new asset, or a new way to build one, meant patching the core codebase. Now you can define your own recipes and asset classes and publish them as a **data channel** on any web host, then share them or subscribe to others', with no changes to core refgenie. Channels are simple collections of YAML files hosted on GitHub Pages or any web server.

See [Use data channels](data_channels.md) and [Set up your own data channel](set_up_data_channel.md).

### Cloud-friendly asset retrieval

Legacy refgenie distributed every asset as a bundled `.tgz` archive, so pulling anything -- including from remote cloud storage -- meant fetching a whole archive to get one file. Asset classes now declare how they are delivered via `serving_modes`:

| Mode | Description |
|---|---|
| `file` | Serve individual files at their own URLs (default) |
| `archive` | Serve a compressed `.tgz` tarball for bulk download |
| `none` | Metadata only -- no data hosting |

With `file` mode you retrieve individual asset files directly from their own URLs, including straight from cloud storage, fetching just the one file you need.

See [Serving modes](serving_modes.md) and [Configure remote storage](remotes.md).

### Upgraded genome identifiers

Legacy refgenie tracked genomes by arbitrary names, so the same assembly could go by different labels and its contents were never verified. Genomes and sequences are now identified by their GA4GH content digest, with human-readable aliases like `hg38` layered on top. Identity is reproducible, and identical sequences are stored only once. The alias system supports both local (RefgetStore-based) and server-side (SQL) backends, so aliases resolve consistently across environments.

See [Manage aliases](alias.md).

### Local dashboard

Managing your local assets was previously CLI-only, with no way to see what you had at a glance. Now `refgenie dash` opens a local web app that lists your genomes, aliases, and assets -- including remote assets from subscribed servers -- with per-genome detail pages, and no server to deploy.

See [Use the dashboard](dash.md).

### AI-assistant access

There was previously no way for an AI tool to see your genome assets. Refgenie now ships a built-in MCP server: point Claude (or another MCP client) at it and ask questions like "what genomes do I have?" or "show the bowtie2 assets for hg38" in plain language. It reads your local database read-only and never modifies anything.

See [Connect AI assistants to your refgenie database](mcp.md).

### Standards-based interoperability

Assets were previously reachable only through refgenie's own API. The built-in server now also exposes standard **GA4GH DRS** endpoints, so any DRS-aware workflow tool can consume refgenie assets directly.

See [Run a server](server/README.md).

### Direct sequence retrieval

Legacy `refgenie seek` returned only a file *path* -- to read the actual sequence you opened the FASTA yourself. Now `refgenie getseq` pulls a sequence or subsequence straight out of the content-addressable store, with no FASTA file on disk required.

See [Genome initialization](genome_tutorial.md) for the full workflow.

## Under the hood

These architectural changes are mostly invisible in day-to-day use but matter if you are migrating, deploying a server, or curious how 1.0 is built.

### Unified package

The legacy ecosystem split functionality across three separate packages:

| Legacy package | Role |
|---|---|
| refgenie | CLI and asset management |
| refgenconf | Configuration library |
| refgenieserver | API server |

Refgenie 1.0 consolidates all of this into one installable package. Server functionality is available with `pip install refgenie[server]`.

### Database-backed configuration

The legacy YAML-based genome configuration file (`genome_config.yaml`) is replaced by a database. By default, refgenie uses a local SQLite database that requires zero configuration. For multi-user or server deployments, PostgreSQL is supported. See [Database backends](database.md) for setup details.

### RefgetStore for sequence data

When you initialize a genome with `refgenie genome init --fasta genome.fa`, refgenie loads all sequences into a local **RefgetStore** -- a content-addressable store where every sequence is identified by its GA4GH digest. This underpins the deduplication, digest-based identity, and `getseq` retrieval described above, and chromosome sizes are derived from sequence collection metadata rather than samtools.

### Built-in server

The `refgenie serve` command starts an API server directly -- no separate refgenieserver package needed. The server exposes REST API endpoints, GA4GH DRS endpoints, and data channel aggregation.

See [Run a server](server/README.md) for setup instructions.

### Staging and remote storage

Assets are **staged** before serving or pushing to cloud storage. Staging prepares assets according to their serving modes -- either as symlinks (file mode, no disk cost) or as `.tgz` tarballs (archive mode). Staged assets can then be pushed to remote storage like S3.

See [Stage assets for serving](staging.md) and [Configure remote storage](remotes.md).

## Migration from legacy

You do not have to rebuild everything. Legacy v0.3/v0.4 installs can be migrated automatically with the [`refgenie-upgrade`](upgrade.md) tool, which reads your legacy `genome_config.yaml`, resolves genome digests to GA4GH format, and registers your existing assets in a fresh refgenie1 database:

```bash
pip install refgenie-upgrade
refgenie-upgrade /path/to/genome_config.yaml
```

See [Upgrading from legacy](upgrade.md) for the full walkthrough.
