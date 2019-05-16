
# <img src="img/refgenie_logo.svg" class="img-header"> genome index manager

[![PEP compatible](http://pepkit.github.io/img/PEP-compatible-green.svg)](http://pepkit.github.io)

## What is refgenie?

Refgenie is full-service reference genome resource manager. It provides a command-line and python interface to download pre-built genome assets like indexes from a central server. It can also build custom genome indexes for any `fasta` file, maintaining a standardized folder structure for reference genome files.

## What makes refgenie better?

Refgenie specifies a *standard* folder structure, so that alignment tools can easily swap from one genome to another. There are other similar projects, but Refgenie has a few advantages:

1. **It provides a command-line interface to download resources**. Think of *Refgenie Server* as `GitHub` for reference genomes. You just type `refgenie pull -g hg38 -a bowtie2...` and you have standard indexes.

2. **It's scripted**. In case you need resources for a custom genome, refgenie provides a `build` function, so users can create their own assembly index packages from any genome source.

3. **It includes a python API**. For tool developers, you can do `cfg = refgenie.RefGenConf(file)` and then use a python object to reference standard locations for any genome asset, *e.g.*, `cfg.hg38.bowtie2_index`.



## Example

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

### Building your own indexes and assets for a reference genome


```console
refgenie build --input hg38.fa.gz
```

