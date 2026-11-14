# What's new in refgenie 1.0

Refgenie 1.0 is a rewrite. One package now does the work of the three 0.x packages (refgenie, refgenconf, and refgenieserver). It installs with `pip install "refgenie>=1.0.0a2"`. You can move an existing 0.x setup to 1.0 without rebuilding your assets; see [Upgrade from legacy](upgrade.md).

## For users

- **Your own recipes and asset types.** Asset classes and recipes come from **data channels**, simple collections of YAML files on any web host. Subscribe to others' channels or publish your own. See [Use data channels](data_channels.md) and [Set up your own data channel](set_up_data_channel.md).
- **Single-file downloads.** Servers can serve each asset file at its own URL, including straight from cloud storage, as well as whole-asset `.tgz` archives. See [Serving modes](serving_modes.md).
- **Standard genome identifiers.** Genomes are identified by GA4GH [sequence collection](https://ga4gh.github.io/seqcol-spec/) digests and sequences by refget digests, with aliases like `hg38` on top. See [Manage aliases](alias.md).
- **Sequence retrieval.** `refgenie getseq` reads a sequence or subsequence straight from the local sequence store, with no FASTA file needed. See [Genome initialization](genome_tutorial.md).
- **A local dashboard.** `refgenie dash` opens a web app that shows your genomes, aliases, and assets, plus assets on your subscribed servers. See [Use the dashboard](dash.md).
- **AI-assistant access.** A built-in, read-only MCP server lets Claude or another MCP client answer questions like "what genomes do I have?" See [Connect AI assistants](mcp.md).
- **GA4GH DRS endpoints.** The server exposes standard DRS endpoints, so any DRS-aware workflow tool can fetch refgenie assets. See [Run a server](server/README.md).

## Under the hood

- **Database instead of a YAML config.** Refgenie keeps its records in SQLite by default, or PostgreSQL for shared and server setups. See [Database backends](database.md).
- **RefgetStore.** `refgenie genome init` loads a genome's sequences into a content-addressable store, where identical sequences are kept once. This store powers `getseq` and the `fasta` asset.
- **Built-in server.** `refgenie serve` (from `pip install "refgenie[server]>=1.0.0a2"`) runs the REST API, DRS endpoints, and data channel aggregation. See [Run a server](server/README.md).
- **Staging and remotes.** Assets are staged for serving, as symlinks or `.tgz` archives, and can be pushed to remote storage like S3. See [Stage assets for serving](staging.md) and [Configure remote storage](remotes.md).
