# Refgenie from within Python

Third-party python tools can rely on our Python object for access to refgenie assets. `refgenie` CLI relies on `Refgenie` object that provides a Python API for accessing and managing genome assets.

## Installing

No additional installation is required if you have `refgenie` installed.

## Quick start

Create a `Refgenie` object, which is the package's main entry point. By default it connects to a local SQLite database. Run `refgenie init` to create one, or point at an existing database configuration file with the `database_config_path` argument (or the `REFGENIE_DB_CONFIG_PATH` environment variable).

Refgenie exposes its functionality through *managers*, accessed as properties on the object: `r.asset`, `r.genome`, `r.alias`, `r.recipe`, `r.asset_class`, `r.servers`, `r.sources`, `r.stage`, and `r.configuration`. The `Refgenie` object itself carries only cross-manager operations, such as `pull`, `getseq`, `populate`, and `build_asset`. There are intentionally no CLI-style delegating wrappers: the `refgenie list ...` command corresponds to `r.asset.table()` in Python, not to an `r.list()` method.

```python
from refgenie import GenomeAlias, Refgenie
r = Refgenie()
```

Now, you can interact with it:

```python
print(r)
```

Use this to show all available remote assets (requires a subscription to a refgenie v4 server):

```python
r.servers.assets_table()
```

In a tool, you're probably most interested in using refgenie to locate reference genome assets. Refgenie groups related files into *asset groups* (such as `fasta` or `bowtie2_index`); within a group, a specific build is an *asset*. Use `r.asset.seek(...)` to get a local file path. For example:

```python
# identify genome (perhaps provided by user)
genome = "hg38"

# get the local path to bowtie2 indexes:
bt2idx = r.asset.seek(r.alias.resolve(GenomeAlias(genome)), "bowtie2_index")

# run bowtie2...
```

This enables you to write python software that will work on any computing environment without having to worry about passing around brittle environment-specific file paths. See [the Refgenie tutorial](notebooks/refgenie.ipynb) for a more comprehensive example of how to work with refgenie as a tool developer.

