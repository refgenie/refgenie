# Installing refgenie

Refgenie needs Python 3.11 or higher. Version 1.0 is a pre-release, so install it from [PyPI](https://pypi.org/project/refgenie/) with `--pre`:

```bash
pip install --pre refgenie
```

Check that it worked:

```bash
refgenie --version
```

## Optional extras

The base install gives you the CLI and the Python API. Add extras for more:

| Extra | Install | What it enables |
| --- | --- | --- |
| `dash` | `pip install --pre "refgenie[dash]"` | `refgenie dash`, the local web UI |
| `server` | `pip install --pre "refgenie[server]"` | `refgenie serve`, the REST API server (includes `dash`) |
| `mcp` | `pip install --pre "refgenie[mcp]"` | `refgenie-mcp`, the MCP server for AI assistants |
| `snakemake` | `pip install --pre "refgenie[snakemake]"` | Building many assets at once with Snakemake |

To install from a source checkout for development, see [Development setup](development.md).

Upgrading from refgenie 0.x? See [Upgrade from legacy](upgrade.md).

## Next steps

[Configure refgenie](configuration.md), or follow the [CLI tutorial](cli_tutorial.md).
