# Use refgenie with an AI assistant

AI assistants such as Claude Code, Claude Desktop and Cursor can work with refgenie in two ways. Use either or both:

1. **The assistant runs refgenie for you.** It runs the CLI, pulls and builds assets, and puts asset paths into your scripts. Point it at the refgenie skill file.
2. **The assistant looks things up in your refgenie database.** It answers questions like "what genomes do I have?" without changing anything. Connect the built-in MCP server.

## Option 1: the skill file

Refgenie publishes an instruction file for AI agents at **<https://docs.refgenie.org/SKILL.md>**. It covers how to `pull`, `build` and `seek`, initialize genomes, retrieve sequences and use the Python API. A second file, <https://api.refgenie.org/SKILL.md>, covers using a refgenie server over HTTP with nothing installed.

Give any assistant that can fetch a URL or run shell commands an instruction like:

> Read https://docs.refgenie.org/SKILL.md, then use refgenie to download the bowtie2 index for hg38 and give me the path.

> Following https://docs.refgenie.org/SKILL.md, initialize my custom genome from `~/data/mygenome.fa.gz` and build its FASTA and bwa indexes.

> Use refgenie (see https://docs.refgenie.org/SKILL.md) to replace the hardcoded genome paths in this alignment script.

!!! tip "Save it in your project"

    If you use refgenie often, save the file into your project so the assistant finds it without being told:

    ```bash
    curl -o SKILL.md https://docs.refgenie.org/SKILL.md
    ```

    If your assistant supports installed skills, you can also put it in that folder, such as a Claude Code skills directory.

!!! note "This option can change your data"

    The assistant runs real refgenie commands, including ones that download or build assets. Review what it proposes, as you would for any command an assistant runs for you. For read-only use, choose the MCP server.

## Option 2: the MCP server

Refgenie includes an **MCP server** that gives an assistant read-only access to your local refgenie database. Install the `mcp` extra and connect it to Claude Code:

```bash
pip install "refgenie[mcp]>=1.0.0a1"
claude mcp add refgenie refgenie-mcp
```

Then ask things like:

> What genomes do I have in refgenie?

> What assets are available for hg38?

> Do I have any mouse genomes?

> How do hg38 and GRCh38 compare?

For Claude Desktop, the full list of tools, and the public server's MCP endpoint, see [Connect AI assistants (MCP)](mcp.md).

## Which to use

| You want to... | Use |
|---|---|
| Have the assistant download, build or wire up assets | The skill file |
| Ask about genomes and assets you already have | The MCP server |
| Keep the assistant read-only | The MCP server |
| Work in a tool that can only fetch URLs | The skill file |

They work well together: the MCP server lets the assistant see your database, and the skill file shows it how to act on it.

## See also

- [The refgenie skill file](https://docs.refgenie.org/SKILL.md)
- [Connect AI assistants (MCP)](mcp.md)
- [CLI tutorial](cli_tutorial.md): the commands the assistant runs, explained for people
