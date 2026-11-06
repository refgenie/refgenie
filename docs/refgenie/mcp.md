# Connect AI assistants (MCP)

Refgenie includes an [MCP](https://modelcontextprotocol.io/) (Model Context Protocol) server that lets AI assistants such as Claude query your refgenie database. Ask "What genomes do I have?" or "Show me the bowtie2 assets for hg38" and the assistant looks up the answer in your refgenie.

The MCP server is read-only. It queries genomes, assets, recipes and sequence collections.

!!! info "Prerequisites"

    - Refgenie installed with the `mcp` extra (`pip install "refgenie[mcp]>=1.0.0a1"`; the `server` extra includes it too)
    - A refgenie database with at least one genome
    - An MCP client, such as Claude Desktop or Claude Code

## Configure your AI client

The MCP server is the `refgenie-mcp` command. It runs as a local process and talks to your AI client over stdio. It needs the `mcp` extra; without it, `refgenie-mcp` exits with an import error.

### Claude Code

```console
claude mcp add refgenie refgenie-mcp
```

Claude Code can now read your refgenie database in every conversation.

### Claude Desktop

Open **Settings > Developer > Edit Config** and add the refgenie server to `mcpServers`:

```json
{
  "mcpServers": {
    "refgenie": {
      "command": "refgenie-mcp"
    }
  }
}
```

If `refgenie-mcp` is not on your PATH, use its full path, which `which refgenie-mcp` prints:

```json
{
  "mcpServers": {
    "refgenie": {
      "command": "/home/you/.local/bin/refgenie-mcp"
    }
  }
}
```

Restart Claude Desktop after saving the config.

### Other MCP clients

Any MCP client that supports stdio can run `refgenie-mcp`. It takes no arguments and uses your default refgenie configuration.

### The public server

The public server offers the same tools over Streamable HTTP at <https://api.refgenie.org/mcp>, for questions about its catalog rather than your own:

```console
claude mcp add --transport http refgenie-public https://api.refgenie.org/mcp
```

Any refgenie server you run with `refgenie serve` has the same `/mcp` endpoint.

## Available tools

All tools are read-only. They give the same answers as the server's HTTP API.

| Tool | Description |
|------|-------------|
| `catalog_summary` | Count genomes, asset groups and assets, and list the species that have assets |
| `list_genomes` | List genomes, most assets first; filter by `has_assets` and `taxon_id` |
| `search_genomes` | Search genomes by digest, alias, any species name (scientific, common or other common names, e.g. brown rat), other genome names (GRCh38), description or accession (substring matching); same filters |
| `get_genome` | Get one genome by `alias` or `digest`, with its FHR metadata and a `species` object holding all its species' names, other common names included |
| `list_assets` | List assets, optionally of one genome (`alias` or `digest`) and/or one `asset_class` |
| `get_asset` | Get one asset by digest: names, seek keys, parents, children, and how to get its files |
| `list_asset_classes` | List asset classes with their seek key names |
| `list_recipes` | List recipes, optionally those that build one `output_asset_class` |
| `get_recipe` | Get one recipe in full: inputs, command templates, docker image, and the build command |
| `lookup_digest` | Say what a digest is: a genome, an asset, or nothing |
| `get_genome_metadata` | Get sequence collection metadata (number of sequences, total length, source) |
| `get_genome_sequences` | Get sequence-level data: names, lengths, and sequence digests |
| `compare_genomes` | Compare two genomes using seqcol comparison |

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

Over stdio, `get_asset` returns `local_path`, where the files are on this machine. On a server (`/mcp`), it returns `downloads`, the archive and file-list URLs relative to the server. Both include `pull_command`, the `refgenie pull` line for the asset.

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

The assistant calls `get_asset` with the digest and returns the full asset record, including its local path, seek keys, any parent/child relationships, and the command to pull it.

**Searching across genomes:**

> "Do I have any mouse genomes?"

The assistant calls `search_genomes` with `query="mouse"` and returns any genomes matching that term in any species name, alias, or description.

**Comparing genomes:**

> "How do hg38 and GRCh38 compare?"

The assistant calls `compare_genomes` with both identifiers and returns the seqcol comparison showing which sequences match, differ, or are unique to each.

!!! success "Key points"

    - The `refgenie-mcp` command runs a local MCP server over stdio
    - All tools are read-only queries against your refgenie database
    - Name a genome by `alias` or by `digest`
    - List tools return 20 small rows per page; follow `next_page_token` for more
    - `refgenie-mcp` takes no arguments and uses your default refgenie configuration
    - Servers offer the same tools at `/mcp`, including <https://api.refgenie.org/mcp>
