
# <img src="img/refgenie_logo.svg" class="img-header"> reference genome manager

[![PEP compatible](https://pepkit.github.io/img/PEP-compatible-green.svg)](https://pepkit.github.io)
[![PyPi](https://img.shields.io/pypi/v/refgenie.svg)](https://pypi.org/project/refgenie/)

## What is refgenie?

Refgenie manages storage, access, and transfer of reference genome resources. It provides command-line and Python interfaces to *download* pre-built reference genome "assets", like indexes used by bioinformatics tools. It can also *build* assets for custom genome assemblies. Refgenie provides programmatic access to a standard genome folder structure, so software can swap from one genome to another.

## What makes refgenie better?

1. **It provides a command-line interface to download individual resources**. Think of it as `GitHub` for reference genomes. You just type `refgenie pull hg38/bwa_index`.

2. **It's scripted**. In case you need resources *not* on the server, such as for a custom genome, you can `build` your own: `refgenie build custom_genome/bowtie2_index`.

3. **It simplifies finding local asset locations**. When you need a path to an asset, you can `seek` it, making your pipelines portable across computing environments: `refgenie seek hg38/salmon_index`.

4. **It includes a python API**. For tool developers, you use `rgc = refgenconf.RefGenConf("genomes.yaml")` to get a Python object with paths to any genome asset, *e.g.*, `rgc.seek("hg38", "kallisto_index")`.


## Quick example

### Install and initialize

Refgenie keeps track of what's available using a configuration file initialized by `refgenie init`:

```console
pip install --user refgenie
export REFGENIE='genome_config.yaml'
refgenie init -c $REFGENIE
```

### Download indexes and assets for a remote reference genome

Use `refgenie pull` to download pre-built assets from a remote server. View available remote assets with `listr`:

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

Refgenie assets are scripted, so if what you need is not available remotely, you can use `build` it locally:


```console
refgenie build mygenome/bwa_index --fasta mygenome.fa.gz
```

See [further reading on building assets](build.md).

### Retrieve paths to refgenie-managed assets

Once you've populated your refgenie with a few assets, use `seek` to retrieve their local file paths:

```console
refgenie seek mm10/bowtie2_index
```

This will return the path to the particular asset of interest, regardless of your computing environment. This gives you an ultra-portable asset manager! See [further reading on retrieving asset paths](seek.md).

If you want to read more about the motivation behind refgenie and the software engineering that makes refgenie work, proceed next to the [overview](overview.md).
