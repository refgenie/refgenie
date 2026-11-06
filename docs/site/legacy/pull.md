# Download pre-built reference genome assets

## Introduction

Use the `refgenie` command-line interface to download and organize genome assets. You do this by simply running `refgenie` from the command line.

The `listr` command *lists remote assets* to see what's available:

```console
refgenie listr
```

The `pull` *downloads* the specific asset of your choice:

```console
refgenie pull GENOME/ASSET
```
Where `GENOME` refers to a genome key (*e.g.* hg38) and `ASSET` refers to one or more specific asset keys (*e.g.* bowtie2_index). For example:

```console
refgenie pull hg38/bowtie2_index
```

You can also pull many assets at once:

```console
refgenie pull --genome mm10 bowtie2_index hisat2_index
```


To see more details, consult the usage docs by running `refgenie pull --help`.

That's it! Easy.

## Downloading manually

You can also browse and download pre-built `refgenie` assemblies manually at [refgenomes.databio.org](http://refgenomes.databio.org).

<!--
## Older builds (deprecated)

For earlier versions of refgenie, there was no `refgenie pull` command, so users would just download the pre-indexed archives and use those. These archives are built for common genomes. These are `tar gzipped` files, so **you will need to unarchive them after downloading**. The complete collection is listed at [http://big.databio.org/refgenomes/](http://big.databio.org/refgenomes/):

## Mirror 1:

* Human: [hg38.tgz](http://big.databio.org/refgenomes/hg38.tgz), [hg19.tgz](http://big.databio.org/refgenomes/hg19.tgz)
* Rat: [rn6.tgz](http://big.databio.org/refgenomes/rn6.tgz)
* Mouse: [mm10](http://big.databio.org/refgenomes/mm10.tgz)
* Prealignment 'decoy' references: available at [big.databio.org/refgenomes/](http://big.databio.org/refgenomes/) for sequences from [ref_decoy](https://github.com/databio/ref_decoy)

## Mirror 2 (use if mirror 1 is down):

* Human: [hg38.tgz](http://cloud.databio.org/refgenomes/hg38.tgz), [hg19.tgz](http://cloud.databio.org/refgenomes/hg19.tgz)
* Rat: [rn6.tgz](http://cloud.databio.org/refgenomes/rn6.tgz)
* Mouse: mm9 (pending), [mm10](http://cloud.databio.org/refgenomes/mm10.tgz)


## Example to download and unarchive hg38_chr22

```
wget http://big.databio.org/refgenomes/hg38_chr22.tgz
tar -xf hg38_chr22.tgz
``` -->
