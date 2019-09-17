# Asset naming scheme

Each asset is defined by *four namespace components*:
 
1. genome name
2. asset name 
3. tag name
4. seek key

The execution of all `refgenie` commands will require at least the first one (genome name) to be specified, but most will also require the asset name. 


## Asset registry paths

The most convenient way to provide this information on the command line is *asset registry path*: `genome/asset.seek_key:tag`, e.g. `hg38/fasta.fai:default`. Yes, that's a lot of typing if you want to be explicit, but `refgenie` makes usage of asset registry paths easy with a system of defaults, such that all the commands below return the same path:

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

## Arguments

Alternatively, you can specify all of these namespace components as command line arguments:

```console
refgenie seek -g rCRSd -a fasta -t default 
```

