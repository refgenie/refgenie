
# <img src="img/refgenie_logo.svg" class="img-header"> reference genome manager

[![PEP compatible](https://pepkit.github.io/img/PEP-compatible-green.svg)](https://pepkit.github.io)


## What is refgenie?

Refgenie is full-service reference genome manager that organizes storage, access, and transfer of reference genomes. It provides command-line and python interfaces to download pre-built reference genome "assets" like indexes used by bioinformatics tools. It can also build assets for custom genome assemblies. Refgenie provides programmatic access to a standard genome folder structure, so software can swap from one genome to another.

## What makes refgenie better?

1. **It provides a command-line interface to download individual resources**. Think of it as `GitHub` for reference genomes. You just type `refgenie pull -g hg38 -a bwa_index`.

2. **It's scripted**. In case you need resources *not* on the server, such as for a custom genome, you can `build` your own: `refgenie build -g custom_genome -a bowtie2_index`.

3. **It simplifies finding local asset locations**. When you need a path to an asset, you can `seek` it, making your pipelines portable across computing environments: `refgenie seek -g hg38 -a salmon_index`.

4. **It includes a python API**. For tool developers, you use `cfg = refgenie.RefGenConf("genomes.yaml")` to get a python object with paths to any genome asset, *e.g.*, `cfg.get_asset("hg38", "kallisto_index")`.


## Quick example

### Install and initialize

```console
pip install --user refgenie
export REFGENIE='genome_config.yaml'
refgenie init -c $REFGENIE
```

### Download indexes and assets for a remote reference genome

First, view available remote assets:

```console
refgenie listr
```

Response:
```console
Querying available assets from server: http://refgenomes.databio.org/v2/assets
Remote genomes: mouse_chrM2x, rCRSd
Remote assets:
        mouse_chrM2x/   bowtie2_index:default, fasta.chrom_sizes:default, fasta.fai:default, fasta:default
               rCRSd/   bowtie2_index:default, fasta.chrom_sizes:default, fasta.chrom_sizes:test, fasta.fai:default, fasta.fai:test, fasta:default, fasta:test
```

Next, pull one:

```console
refgenie pull rCRSd/bowtie2_index
```

Response:
```console
'rCRSd/bowtie2_index:default' archive size: 116.8KB
Downloading URL: http://staging.refgenomes.databio.org/v2/asset/rCRSd/bowtie2_index/archive ... 
```

See [further reading on downloading assets](pull.md).

### Build your own indexes and assets for a custom reference genome


```console
refgenie build mygenome/bwa_index --fasta mygenome.fa.gz
```

See [further reading on building assets](build.md).

### Retrieve paths to refgenie-managed assets

Once you've populated your refgenie with a few assets, it's easy to get paths to them:

```console
refgenie seek mm10/bowtie2_index
```

This will return the path to the particular asset of interest, regardless of your computing environment. This gives you an ultra-portable asset manager! See [further reading on retrieving asset paths](seek.md).

If you want to read more about the motivation behind refgenie and the software engineering that makes refgenie work, proceed next to the [overview](overview.md).
