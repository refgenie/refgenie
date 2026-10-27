# How to connect AI assistants to your refgenie database

Refgenie includes a built-in [MCP](https://modelcontextprotocol.io/) (Model Context Protocol) server that lets AI assistants like Claude query your local refgenie database. Once configured, you can ask your AI assistant things like "What genomes do I have?" or "Show me the bowtie2 assets for hg38" and it will look up the answers directly from your refgenie instance.

The MCP server is read-only -- it can query genomes, assets, recipes, and sequence collections, but it never modifies your data.

!!! info "Prerequisites"

    - You have refgenie installed (`pip install refgenie`)
    - You have a refgenie database initialized with at least one genome
    - You have an MCP-compatible AI client (Claude Desktop, Claude Code, etc.)

## Configure your AI client

The MCP server runs as a local process using stdio transport. You point your AI client at the `refgenie-mcp` command, and it handles the rest.

### Claude Code

Add the refgenie MCP server with a single command:

```console
claude mcp add refgenie refgenie-mcp
```

That's it. Claude Code will now have access to your refgenie database in every conversation.

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

If `refgenie-mcp` is not on your system PATH, use the full path to the executable. You can find it with:

```console
which refgenie-mcp
```

Then use that full path in the config:

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

Any MCP-compatible client that supports stdio transport can use the refgenie server. The command is `refgenie-mcp` with no arguments. It reads from stdin and writes to stdout using the MCP protocol.

## Available tools

The MCP server exposes these read-only tools:

| Tool | Description |
|------|-------------|
| `list_genomes` | List all genomes with aliases, species, and description |
| `list_asset_classes` | List registered asset classes with seek keys and serving modes |
| `list_recipes` | List registered recipes with output asset class and inputs |
| `list_assets` | List assets, optionally filtered by genome and/or asset class |
| `search_genomes` | Search genomes by species name, alias, or description (substring matching) |
| `get_genome` | Get detailed info for a genome by alias or digest, including its asset groups |
| `get_asset` | Get detailed asset info by digest, including seek keys, parents, and children |
| `lookup_digest` | Universal digest lookup -- tries genome first, then asset |
| `get_genome_metadata` | Get sequence collection metadata (number of sequences, total length, source) |
| `get_genome_sequences` | Get sequence-level data: names, lengths, and sequence digests |
| `compare_genomes` | Compare two genomes using seqcol comparison |

Genomes can be referenced by alias (like `hg38`) or by digest. The server resolves aliases automatically.

## Example interactions

Once the MCP server is configured, you can ask your AI assistant questions in natural language. Here are a few examples of what you might ask and what happens behind the scenes.

**Browsing your genomes:**

> "What genomes do I have in refgenie?"

The assistant calls `list_genomes` and returns a summary of all your genomes with their aliases and species.

**Finding assets for a genome:**

> "What assets are available for hg38?"

The assistant calls `list_assets` with `genome="hg38"` and lists all assets built or downloaded for that genome.

**Looking up a specific asset:**

> "Get me the details on this asset: a1b2c3d4e5"

The assistant calls `get_asset` with the digest and returns the full asset record, including its path, seek keys, and any parent/child relationships.

**Searching across genomes:**

> "Do I have any mouse genomes?"

The assistant calls `search_genomes` with `query="mouse"` and returns any genomes matching that term in their species name, alias, or description.

**Comparing genomes:**

> "How do hg38 and GRCh38 compare?"

The assistant calls `compare_genomes` with both identifiers and returns the seqcol comparison showing which sequences match, differ, or are unique to each.

!!! success "Key points"

    - The `refgenie-mcp` command runs a local MCP server over stdio
    - All tools are read-only queries against your refgenie database
    - Genomes can be referenced by alias or digest -- the server resolves both
    - No arguments or configuration flags are needed; `refgenie-mcp` picks up your default refgenie configuration
