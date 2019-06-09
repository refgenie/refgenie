
# <img src="img/refgenie_logo.svg" class="img-header"> reference genome manager

[![PEP compatible](http://pepkit.github.io/img/PEP-compatible-green.svg)](http://pepkit.github.io)


## What is refgenie?

Refgenie is full-service reference genome manager. It provides command-line and python interfaces to download pre-built genome assets like indexes. It can also build assets for custom genomes.

## What makes refgenie better?

Refgenie provides programmatic access to a standard genome folder structure, so that software can easily swap from one genome to another. There are other similar projects, but Refgenie has a few advantages:

1. **It provides a command-line interface to download individual resources**. Think of it as `GitHub` for reference genomes. You just type `refgenie pull -g hg38 -a bowtie2_index`.

2. **It's scripted**. In case you need resources *not* on the server, such as for a custom genome, refgenie provides a `build` function to create your own: `refgenie build -i custom.fa.gz -a bowtie2_index`.

3. **It includes a python API**. For tool developers, you use `cfg = refgenie.RefGenConf("genomes.yaml")` to get a python object with paths to any genome asset, *e.g.*, `cfg.hg38.bowtie2_index`.



## Quick example

### Downloading indexes and assets for a reference genome


```console
refgenie pull --genome hg38 --asset bowtie2_index
```

Response:
```console
Pulling assets 'bowtie2_index' from genome 'hg38'.
Downloading URL: http://.../asset/hg38/bowtie2_index
Download complete: /path/to/genome/hg38/bowtie2_index.tgz
Unpacked archive at: /path/to/genome/hg38/bowtie2_index
```

Pull many assets at once:
```console
refgenie pull --genome mm10 --asset kallisto TSS_enrichment mappability
```

See [further reading on downloading assets](download.md).

### Building your own indexes and assets for a reference genome


```console
refgenie build --input hg38.fa.gz --asset bowtie2_index
```

See [further reading on building assets](build.md).

If you want to read more about the motivation behind refgenie and the software engineering that makes refgenie work, proceed next to the [overview](overview.md).