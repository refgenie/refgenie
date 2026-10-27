# Installing refgenie

Install from [PyPI](https://pypi.org/project/refgenie/). Refgenie 1.0 is currently a pre-release, so tell `pip` to accept it:

```bash
pip install --pre refgenie
```

## Verify installation

```bash
refgenie --version
```

## Optional extras

Install with dashboard support (local web UI for browsing assets):

```bash
pip install --pre "refgenie[dash]"
```

Install with server support (includes dash plus background scheduling):

```bash
pip install --pre "refgenie[server]"
```

### Development install

To install from a local clone for development:

```bash
pip install -e ".[dash]"
```

## Upgrading from legacy refgenie

If you are migrating from legacy refgenie (v0.3/v0.4), you do not need to rebuild your assets. After installing, use the `refgenie-upgrade` tool to migrate your existing data and configuration:

```bash
pip install refgenie-upgrade
refgenie-upgrade /path/to/genome_config.yaml
```

See the [upgrade tutorial](upgrade.md) for the full walkthrough.

## Next steps

Configure refgenie!

For more information on how to configure refgenie, see the [Configuration documentation](configuration.md).