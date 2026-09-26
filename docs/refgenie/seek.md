# Retrieve paths to *local* assets

Once you've assembled a few assets, either by downloading or by building them, you'll be able to use `refgenie seek` to retrieve the paths. It's quite simple, really -- say you've built the `bowtie2_index` and `fasta` assets for `hg38`. If you type:

```console
refgenie seek hg38/bowtie2_index
```

You'll get back the absolute path on your system to the `bowtie2_index` asset, something like:

```console
/path/to/genomes_folder/hg38/bowtie2_index
```

## Seek keys

Refgenie `seek` uses a concept called *seek keys*. Let's formally define how these differ from an asset:

* **asset**: a folder consisting of one or more files related to a specific reference genome
* **seek key**: a string identifier for a particular file or folder contained within an asset that can be located with `refgenie seek`.

An asset can have multiple seek keys. When an asset recipe is defined, the author records which files in the asset are potentially interesting for a user to retrieve a path to. These files are added to the recipe as seek keys.

Many assets have only a single seek key, because the asset really only consists of a single thing. For example, a bowtie2_index asset has only 1 thing you'd need to *seek*: the folder containing the indexes. But other assets can have multiple seek keys. For example, the `fasta` asset provides not only a fasta file, but a `.fai` fasta index as well as a *chrom_sizes* file that is required by some tools (notably, UCSC tools). These files are really tightly coupled with a particular fasta file and it doesn't make sense to make separate assets for each of them, and this is where seek keys come in.

Once you have built a `fasta` asset, you'll be able to retrieve not just the fasta file, but also the `fai`, and `chrom_sizes` files.  These will show up when you type `refgenie list` as `fasta.fai` and `fasta.chrom_sizes` seek keys.

## Seek keys make things portable

Why is this useful? Instead of writing your shell script to take a path to the `chrom_sizes` file and the `bowtie2_index` path, you just take the genome name, and use `refgenie seek` to find the correct path. Not only have you simplified the interface to your analysis, but both the script itself *and your call to it* are now portable.

In other words, you just write:

```console
python script.py --genome hg38
```

and then use refgenie to *seek* any assets you need. Previously, you would have done:

```console
python script.py --chrom_sizes /local/path/to/hg38/chrom_size \
  --bt2_index /local/path/to/hg38/bowtie2_index
```

This second command obviously differs by computing environment, whereas the first does not.

## Example

```console
refgenie seek hg38/bowtie2_index
refgenie seek hg38/fasta.chrom_sizes
refgenie seek hg38/fasta
```
