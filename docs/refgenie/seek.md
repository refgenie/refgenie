# Retrieve paths to local assets

`refgenie seek` prints the path to an asset you have pulled or built:

```console
$ refgenie seek hg38/bowtie2_index
/path/to/genomes/alias/hg38/bowtie2_index/2.5.4/hg38
```

## Seek keys

An *asset* is a folder of one or more files for one genome. A *seek key* names one file or folder inside it. The asset class decides which seek keys an asset has.

Many assets have one seek key. A `bowtie2_index` has just the index prefix. Others have several: a `fasta` asset holds the FASTA file, its `.fai` index, and a `chrom_sizes` file, which belong together, so they are one asset with three seek keys. Add the seek key after a period:

```console
refgenie seek hg38/fasta              # the FASTA file
refgenie seek hg38/fasta.fai          # the FASTA index
refgenie seek hg38/fasta.chrom_sizes  # the chromosome sizes file
```

`refgenie list` shows each asset's seek keys. See [Asset registry paths](asset_registry_paths.md) for the full path syntax and its defaults.

## Why this helps

A script that takes a genome name and seeks its files runs the same on any computer:

```console
python script.py --genome hg38
```

A script that takes file paths needs different arguments on each computer:

```console
python script.py --chrom_sizes /local/path/to/hg38.chrom.sizes \
  --bt2_index /local/path/to/hg38/bowtie2_index
```

Add `-e` (`--check-exists`) to fail when the file is missing. To get the address of an asset on a server instead of a local path, use `refgenie seekr`.
