# Asset registry paths

Each asset is defined by four components:

1. genome name
2. asset name
3. tag name
4. seek key

All `refgenie` commands require a genome name, and most also require an asset name. Tag and seek keys are used only when needed and otherwise use sensible defaults.

The most convenient way to provide this information on the command line is with an *asset registry path*, which take this form:

```console
genome/asset.seek_key:tag
```

For example, `hg38/fasta.fai:default`. Yes, that's a lot of typing if you want to be explicit, but `refgenie` makes usage of asset registry paths easy with a system of defaults, such that all the commands below return the same path:

```console
$ refgenie seek rCRSd/fasta
path/to/genomes/archive/rCRSd/fasta/default/rCRSd.fa

$ refgenie seek rCRSd/fasta.fasta
path/to/genomes/archive/rCRSd/fasta/default/rCRSd.fa

$ refgenie seek rCRSd/fasta.fasta:default
path/to/genomes/archive/rCRSd/fasta/default/rCRSd.fa
```

How did it work?

- **default tag** is determined by `default_tag` pointer in the config
- **seek_key** defaults to the name of the asset

## Using arguments instead of registry paths

Alternatively, you can specify all of these namespace components as command line arguments:

```console
refgenie seek -g rCRSd -a fasta -t default
```

One advantage of this approach is that it allows you to refer to multiple assets belonging to the same genome.
