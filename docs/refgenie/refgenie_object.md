# Refgenie from within Python

Python tools can use the `Refgenie` object to find and manage assets. It comes with the `refgenie` package; nothing else needs installing.

```python
from refgenie import GenomeAlias, Refgenie

r = Refgenie()
```

By default this opens the local SQLite database that `refgenie init` creates. To use another database, pass `database_config_path` or set `REFGENIE_DB_CONFIG_PATH`.

The object groups its features into *managers*: `r.asset`, `r.genome`, `r.alias`, `r.recipe`, `r.asset_class`, `r.servers`, `r.sources` (data channels), `r.stage`, `r.remote`, `r.transfer` (pull), `r.build`, `r.sequence`, and `r.configuration`. Most CLI commands map to one manager method. For example, `refgenie list` is `r.asset.table()`, and `refgenie pull` is `r.transfer.pull(...)`.

Show the assets on your subscribed servers:

```python
r.servers.assets_table()
```

In a tool, you will mostly want local file paths. Refgenie groups related files into *asset groups* (such as `fasta` or `bowtie2_index`); within a group, each build is an *asset*. `r.asset.seek(...)` returns a path:

```python
genome = "hg38"  # perhaps provided by the user

bt2idx = r.asset.seek(r.alias.resolve(GenomeAlias(genome)), "bowtie2_index")
```

Your tool then takes a genome name instead of file paths, so it runs the same on any computer. See [Use refgenie in your pipeline](code_snippets.md) for a function that also pulls a missing asset, and [the Refgenie tutorial](notebooks/refgenie.ipynb) for a longer example.
