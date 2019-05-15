# <img src="img/refgenie_logo.svg" class="img-header"> genome index manager

[![PEP compatible](http://pepkit.github.io/img/PEP-compatible-green.svg)](http://pepkit.github.io)

## What is refgenie?

Refgenie creates a standardized folder structure for reference genome files and indexes. You can download pre-built genomes or build your own for any fasta file.

## What makes refgenie better?

Refgenie provides a **standard folder structure** for reference genome indexes, so that alignment tools can easily swap from one genome to another. Most importantly, Refgenie is **scripted** so that users can create their own assembly index packages from whatever genome source they like.

## Example

Download aligner indexes and other resources for your genome of interest right from the command-line:


```console
refgenie pull --genome hg38 --asset bowtie2
```

Response:
```console
Pulling... Genome: hg38; assets: bowtie2
Downloading... URL: http://.../asset/hg38/bowtie2
Download complete.
Saved as: /path/to/genome/hg38/bowtie2.tar
Unarchived result at: /path/to/genome/hg38/bowtie2_index
```

Pull many assets at once:
```console
refgenie pull --genome mm10 --asset kallisto TSS_enrichment mappability
```

Or, build your own indexes:

```console
refgenie build --input hg38.fa.gz
```