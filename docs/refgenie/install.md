# Installing refgenie

Refgenie needs Python 3.11 or higher. Version 1.0 is a pre-release, so install it from [PyPI](https://pypi.org/project/refgenie/) by version:

```bash
pip install "refgenie>=1.0.0a2"
```

Without the version, `pip install refgenie` gives you the old 0.13, and so does Bioconda. Version 0.13 uses different commands, so the examples in these docs fail with it.

Check that you have 1.0:

```bash
refgenie --version   # refgenie 1.0.0a2
```

Used refgenie 0.x, or copied a command from older docs? Refgenie 1.0 ignores `$REFGENIE` (it reads `REFGENIE_HOME_PATH`), has no `refgenie init -c`, and names assets as paths (`refgenie seek hg38/fasta`), not with `-g hg38 -a fasta`. See [old commands and their 1.0 forms](migration-from-legacy.md#commands) and [Upgrade from 0.x](upgrade.md).

## Optional extras

The base install gives you the CLI and the Python API. Add extras for more:

| Extra | Install | What it enables |
| --- | --- | --- |
| `dash` | `pip install "refgenie[dash]>=1.0.0a2"` | `refgenie dash`, the local web UI |
| `server` | `pip install "refgenie[server]>=1.0.0a2"` | `refgenie serve`, the REST API server (includes `dash` and `mcp`) |
| `mcp` | `pip install "refgenie[mcp]>=1.0.0a2"` | `refgenie-mcp`, the MCP server for AI assistants |
| `snakemake` | `pip install "refgenie[snakemake]>=1.0.0a2"` | Building many assets at once with Snakemake |

To install from a source checkout for development, see [Development setup](development.md).

## Next steps

[Configure refgenie](configuration.md), or follow the [CLI tutorial](cli_tutorial.md).
