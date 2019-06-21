
# <img src="img/refgenie_logo.svg" class="img-header"> reference genome manager

[![PEP compatible](http://pepkit.github.io/img/PEP-compatible-green.svg)](http://pepkit.github.io)


## What is refgenie?

Refgenie is full-service reference genome manager. It provides command-line and python interfaces to download pre-built reference genome "assets" like indexes used by different bioinformatics tools. It can also build assets for custom genome assemblies, and it facilitates systematic organization of, and access to, local genome "assets."

## What makes refgenie better?

Refgenie provides programmatic access to a standard genome folder structure, so that software can easily swap from one genome to another. There are other similar projects, but Refgenie has a few advantages:

1. **It provides a command-line interface to download individual resources**. Think of it as `GitHub` for reference genomes. You just type `refgenie pull -g hg38 -a kallisto_index`.

2. **It's scripted**. In case you need resources *not* on the server, such as for a custom genome, refgenie provides a `build` function to create your own: `refgenie build -i custom.fa.gz -a bowtie2_index`.

3. **It includes a python API**. For tool developers, you use `cfg = refgenie.RefGenConf("genomes.yaml")` to get a python object with paths to any genome asset, *e.g.*, `cfg.hg38.kallisto_index`.

4. **It maintains a repository of local asset locations**. Refgenie maintains a local configuration file with the metadata for each resource that's been downloaded or built.

## Quick example

### Install and initialize

```console
pip install --user refgenie
export REFGENIE='genome_config.yaml'
refgenie init -c $REFGENIE
refgenie listr
```

### Downloading indexes and assets for a reference genome

```console
refgenie pull --genome hg38 --asset bowtie2
```

Response:
```console
Starting pull for 'hg38/bowtie2'
'hg38/bowtie2' archive size: 3.5GB
Downloading URL: http://refgenomes.databio.org/asset/hg38/bowtie2/archive
...
```

Pull many assets at once:
```console
refgenie pull --genome mm10 --asset kallisto TSS_enrichment mappability
```

See [further reading on downloading assets](download.md).

### Building your own indexes and assets for a reference genome


```console
refgenie build --genome hg38 --asset kallisto_index --fasta hg38.fa.gz
```

See [further reading on building assets](build.md).

If you want to read more about the motivation behind refgenie and the software engineering that makes refgenie work, proceed next to the [overview](overview.md).
