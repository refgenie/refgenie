# Installing refgenie

Refgenie needs Python 3.11 or higher. Version 1.0 is a pre-release, so install it from [PyPI](https://pypi.org/project/refgenie/) by version:

```bash
pip install "refgenie>=1.0.0a1"
```

Check that it worked:

```bash
refgenie --version
```

## Optional extras

The base install gives you the CLI and the Python API. Add extras for more:

| Extra | Install | What it enables |
| --- | --- | --- |
| `dash` | `pip install "refgenie[dash]>=1.0.0a1"` | `refgenie dash`, the local web UI |
| `server` | `pip install "refgenie[server]>=1.0.0a1"` | `refgenie serve`, the REST API server (includes `dash`) |
| `mcp` | `pip install "refgenie[mcp]>=1.0.0a1"` | `refgenie-mcp`, the MCP server for AI assistants |
| `snakemake` | `pip install "refgenie[snakemake]>=1.0.0a1"` | Building many assets at once with Snakemake |

To install from a source checkout for development, see [Development setup](development.md).

Upgrading from refgenie 0.x? Refgenie 1.x ignores the `$REFGENIE` variable and uses `REFGENIE_HOME_PATH` and `REFGENIE_DB_CONFIG_PATH` instead. See [Upgrade from legacy](upgrade.md).

## Next steps

[Configure refgenie](configuration.md), or follow the [CLI tutorial](cli_tutorial.md).
