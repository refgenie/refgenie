# <img src="img/refgenie_logo.svg" class="img-header"> genome index manager

[![PEP compatible](http://pepkit.github.io/img/PEP-compatible-green.svg)](http://pepkit.github.io)

## What is refgenie?

Refgenie creates a standardized folder structure for reference genome files and indexes. You can download pre-built genomes or build your own for any fasta file.

## What makes refgenie better?

Refgenie provides a **standard folder structure** for reference genome indexes, so that alignment tools can easily swap from one genome to another. Most importantly, Refgenie is **scripted** so that users can create their own assembly index packages from whatever genome source they like.

## Example

### Downloading aligner indexes and assets for a reference genome:


```console
refgenie pull --genome hg38 --asset bowtie2
```

Response:
```console
Pulling assets 'bowtie2' from genome 'hg38'.
Downloading URL: http://.../asset/hg38/bowtie2
Download complete: /path/to/genome/hg38/bowtie2.tar
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

