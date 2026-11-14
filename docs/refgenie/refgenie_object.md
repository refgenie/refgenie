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

Your tool then takes a genome name instead of file paths, so it runs the same on any computer. See [Use refgenie in your pipeline](code_snippets.md) for a function that also pulls a missing asset.

## Worked example

This example uses a throwaway folder, so it never touches `~/.refgenie`. Set the folder before importing refgenie.

```python
import os, tempfile

os.environ["REFGENIE_HOME_PATH"] = tempfile.mkdtemp()
```

Open refgenie and create its database, like `refgenie init`:

```python
from refgenie import GenomeAlias, Refgenie

r = Refgenie(suppress_migrations=True)
r.database.init()
```

Pull the `fasta` and `bowtie2_index` assets of `rCRSd` (a small genome) from the default server, like `refgenie pull rCRSd/fasta rCRSd/bowtie2_index`:

```python
for group in ("fasta", "bowtie2_index"):
    r.transfer.pull(group, GenomeAlias("rCRSd"))
```

Find the files, like `refgenie seek`:

```python
digest = r.alias.resolve(GenomeAlias("rCRSd"))
r.asset.seek(digest, "fasta")  # the .fa
r.asset.seek(digest, "fasta", seek_key_name="chrom_sizes")  # the .chrom.sizes
r.asset.seek(digest, "bowtie2_index")  # the index prefix
```

Register a genome from a FASTA file, like `refgenie genome init`:

```python
from pathlib import Path

fasta = Path(os.environ["REFGENIE_HOME_PATH"]) / "tiny.fa"
fasta.write_text(">chr1\nACGTACGTACGTACGTACGTACGT\n>chr2\nGGGGCCCCAAAATTTT\n")
tiny, created = r.build.initialize_and_build([GenomeAlias("tiny")], fasta_file_path=fasta)
```

Compare two genomes, like `refgenie compare`. The result is the GA4GH seqcol comparison:

```python
r.genome.compare(digest, tiny)["attributes"]  # a_only, b_only, a_and_b per attribute
```

Add and remove an alias, like `refgenie alias set` and `alias remove`:

```python
r.set_genome_alias(GenomeAlias("tiny_copy"), tiny)
r.alias.remove(GenomeAlias("tiny_copy"))
```

To build more assets, see [Build assets](build.md).
