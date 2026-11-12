# Connect AI assistants (MCP)

Refgenie includes an [MCP](https://modelcontextprotocol.io/) (Model Context Protocol) server that lets AI assistants such as Claude query your refgenie database. Ask "What genomes do I have?" or "Show me the bowtie2 assets for hg38" and the assistant looks up the answer in your refgenie.

The MCP server is read-only. It queries genomes, assets, recipes and sequence collections.

!!! info "Prerequisites"

    - Refgenie installed with the `mcp` extra (`pip install "refgenie[mcp]>=1.0.0a1"`; the `server` extra includes it too)
    - `refgenie init` run once. An empty database works; the tools return empty lists until you pull or build something
    - An MCP client, such as Claude Desktop or Claude Code

## Configure your AI client

The MCP server is the `refgenie-mcp` command. It runs as a local process and talks to your AI client over stdio. It needs the `mcp` extra; without it, `refgenie-mcp` exits with an import error.

### Choose the refgenie folder

`refgenie-mcp` reads `REFGENIE_HOME_PATH`, like the CLI (default `~/.refgenie`). Apps started from a desktop launcher do not see your shell's variables, so set it in the server config.

Claude Code:

```console
claude mcp add --env REFGENIE_HOME_PATH=/data/refgenie refgenie -- refgenie-mcp
```

Claude Desktop (**Settings > Developer > Edit Config**) and other clients that use `mcpServers`; restart the app after saving:

```json
{ "mcpServers": { "refgenie": {
    "command": "/full/path/to/refgenie-mcp",
    "env": { "REFGENIE_HOME_PATH": "/data/refgenie" } } } }
```

Codex:

```toml
# ~/.codex/config.toml
[mcp_servers.refgenie]
command = "/full/path/to/refgenie-mcp"
env = { REFGENIE_HOME_PATH = "/data/refgenie" }
```

`which refgenie-mcp` prints the full path. `refgenie-mcp --help` lists its options; on start it prints the database it opened, so check your client's MCP log to confirm the folder.

### The public server

The public server offers the same tools over Streamable HTTP at <https://api.refgenie.org/mcp>, for questions about its catalog rather than your own:

```console
claude mcp add --transport http refgenie-public https://api.refgenie.org/mcp
```

```json
{ "mcpServers": { "refgenie-public": { "type": "http", "url": "https://api.refgenie.org/mcp" } } }
```

Register both servers under different names; they have the same tools.

Any refgenie server you run with `refgenie serve` has the same `/mcp` endpoint.

## Available tools

All tools are read-only. They give the same answers as the server's HTTP API.

| Tool | Description |
|------|-------------|
| `catalog_summary` | Count genomes, asset groups and assets, and list the species that have assets |
| `list_genomes` | List genomes, most assets first; filter by `has_assets` and `taxon_id` |
| `search_genomes` | Search genomes by digest, alias, any species name (scientific, common or other common names, e.g. brown rat), other genome names (GRCh38), description or accession (substring matching; a description must match at a word start); best matches first; same filters |
| `get_genome` | Get one genome by `alias` or `digest`, with its FHR metadata and a `species` object holding all its species' names, other common names included |
| `list_assets` | List assets, optionally of one genome (`alias` or `digest`) and/or one `asset_class` |
| `get_asset` | Get one asset by digest: names, seek keys, parents, children, and how to get its files |
| `list_asset_classes` | List asset classes with their seek key names |
| `list_recipes` | List recipes, optionally those that build one `output_asset_class` |
| `get_recipe` | Get one recipe in full: inputs, command templates, docker image, and the build command |
| `lookup_digest` | Say what a digest is: a genome, an asset, or nothing |
| `get_genome_sequences` | List a genome's sequences (name, length, digest), paged; filter with `name_contains` |
| `get_sequence` | Get the bases of one region (`name`, `start`, `end`; 0-based, half-open; at most 100,000 bases) without downloading the genome |
| `compare_genomes` | Compare two genomes using seqcol comparison, with a one-sentence summary and a map of renamed sequences |

A tool that takes a genome has two parameters, `alias` (a name like `hg38`) and `digest` (32 characters). Give exactly one.

### Paging

List tools return 20 rows per page by default and at most 100 (`page_size`; a larger value is lowered to 100, and the result shows the size used). Each result is a GA4GH token-paged envelope:

```json
{
  "results": [{"digest": "...", "aliases": ["hg38"], "species_name": "Homo sapiens", "common_name": "Human", "asset_count": 12}],
  "pagination": {"page_size": 20, "next_page_token": "eyJ2IjoxLCJvIjoyMCwiaCI6Ii4uLiJ9", "total": 701},
  "hint": null
}
```

To get the next page, call the tool again with the same filters and `token` set to `next_page_token`. On the last page it is `null`. A token from a query with other filters is refused. When `search_genomes` finds nothing, `hint` suggests what to try instead.

### Getting the files

`get_asset` returns `files`, one entry per file (local paths over stdio, URLs on a server), `archive_url` for archive-served assets, and `pull_command`, the `refgenie pull` line for the asset.

### Errors

A failed call comes back as a tool error that names the problem and the next step, for example: "Genome alias 'hg39' not found. Find genomes with search_genomes(query) or list_genomes()." An unknown argument is an error that lists the arguments the tool takes. For example, `list_assets(genome="hg38")` is refused.

## Example interactions

Here are some questions and the tools the assistant calls to answer them.

**Browsing your genomes:**

> "What genomes do I have in refgenie?"

The assistant calls `list_genomes` and returns a summary of your genomes with their aliases and species, following `next_page_token` if there are more than fit on one page.

**Finding assets for a genome:**

> "What assets are available for hg38?"

The assistant calls `list_assets` with `alias="hg38"` and lists the assets built or downloaded for that genome.

**Looking up a specific asset:**

> "Get me the details on this asset: a1b2c3d4e5"

The assistant calls `get_asset` with the digest and returns the full asset record, including its files, seek keys, any parent/child relationships, and the command to pull it.

**Searching across genomes:**

> "Do I have any mouse genomes?"

The assistant calls `search_genomes` with `query="mouse"` and returns any genomes matching that term in any species name, alias, or description.

**Reading a region:**

> "Show me the first 100 bases of BRCA1 in hg38."

The assistant calls `get_sequence` with `alias="hg38"`, `name="chr17"`, `start=43044294`, `end=43044394`.

**Comparing genomes:**

> "How do hg38 and GRCh38 compare?"

The assistant calls `compare_genomes` with both identifiers and returns the seqcol comparison showing which sequences match, differ, or are unique to each.

!!! success "Key points"

    - The `refgenie-mcp` command runs a local MCP server over stdio
    - All tools are read-only queries against your refgenie database
    - Name a genome by `alias` or by `digest`
    - List tools return 20 small rows per page; follow `next_page_token` for more
    - `refgenie-mcp` reads `REFGENIE_HOME_PATH`
    - Servers offer the same tools at `/mcp`, including <https://api.refgenie.org/mcp>
