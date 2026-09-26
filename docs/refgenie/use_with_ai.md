# How to use refgenie with an AI assistant

AI coding assistants (Claude Code, Claude Desktop, Cursor, and similar tools)
work well with refgenie. There are two complementary ways to combine them, and
you can use either or both:

1. **Let the assistant *operate* refgenie for you** — run the CLI, pull and
   build assets, and wire asset paths into your scripts. Point it at the
   refgenie skill file.
2. **Let the assistant *query* your refgenie database** — answer questions like
   "what genomes do I have?" without touching anything. Use the built-in MCP
   server.

## Option 1: Point your assistant at the skill file

Refgenie publishes a single self-contained instruction file for AI agents at:

**<https://docs.refgenie.org/SKILL.md>**

It teaches an assistant the refgenie commands and concepts it needs — how to
`pull`, `build`, `seek`, initialize genomes, retrieve sequences, and use the
Python API — in one page. You do not need to explain refgenie yourself; you
just direct the assistant to read it.

### Direct the assistant in plain language

With any assistant that can fetch a URL or run shell commands, give an
instruction like:

> Read https://docs.refgenie.org/SKILL.md, then use refgenie to download the bowtie2
> index for hg38 and give me the path.

> Following https://docs.refgenie.org/SKILL.md, initialize my custom genome from
> `~/data/mygenome.fa.gz` and build its FASTA and bwa indexes.

> Use refgenie (see https://docs.refgenie.org/SKILL.md) to add the correct genome
> asset paths to this alignment script instead of the hardcoded ones.

The assistant reads the page, learns the commands, and runs them for you.

!!! tip "Save it once for a project"

    If you work with refgenie often, save the file into your project so the
    assistant picks it up automatically:

    ```bash
    curl -o SKILL.md https://docs.refgenie.org/SKILL.md
    ```

    Coding assistants that read repository context will find and follow it
    without being told each time. You can also drop it into a tool-specific
    location — for example a Claude Code skills directory — if your assistant
    supports installed skills.

!!! note "This mode can change your data"

    Following the skill lets the assistant run real refgenie commands,
    including ones that download or build assets. Review what it proposes, the
    same as you would with any command an assistant runs on your behalf. For a
    read-only experience, use the MCP server below.

## Option 2: Let the assistant query your database (MCP)

Refgenie ships a built-in **MCP server** that gives an AI assistant read-only
access to your local refgenie database. Once connected, you can ask questions
in plain language and the assistant answers from your actual data — it can look
things up, but it never modifies anything.

Connect it in one command with Claude Code:

```bash
claude mcp add refgenie refgenie-mcp
```

Then ask things like:

> What genomes do I have in refgenie?

> What assets are available for hg38?

> Do I have any mouse genomes?

> How do hg38 and GRCh38 compare?

For Claude Desktop configuration, the full list of available tools, and more
examples, see [How to connect AI assistants to your refgenie database](mcp.md).

## Which should I use?

| You want to... | Use |
|---|---|
| Have the assistant download, build, or wire up assets | The skill file (Option 1) |
| Ask questions about genomes and assets you already have | The MCP server (Option 2) |
| Keep the assistant strictly read-only | The MCP server (Option 2) |
| Work in an editor/agent that only fetches URLs | The skill file (Option 1) |

The two combine naturally: connect the MCP server so the assistant can *see*
your database, and hand it the skill file so it knows how to *act* on it.

## See also

- [The refgenie skill file](https://docs.refgenie.org/SKILL.md) — the instruction page for AI agents
- [Connect AI assistants (MCP)](mcp.md) — full MCP setup and tool reference
- [CLI tutorial](cli_tutorial.md) — the commands the assistant runs, explained for humans
