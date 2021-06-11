# Genomes compatibility

## Motivation

Many genomic analyses require less stringent comparison: simple compatibility between assets that are not necessarily identical.
For example, if we wanted to annotate genomic regions based on aligned BAM file using feature annotations we would need the reads to share just the coordinate system (the sequence lengths and names). Importantly, it does not require sequence identity at all. In this case, we would like to confirm that the given `bowtie2_index` asset is at least compatible with the coordinate system of feature annotation file.

To sum up, we need a more detailed comparison that may ignore sequences but allow to compare other aspects of the reference genome.

## Solution

Refgenie facilitates such fine-grained comparison of genomes via `refgenie compare` command. A useful "byproduct" of genome initialization described in [Use aliases how-to guide](aliases.md) is a JSON file with annotated sequence digests, that are required to compare FASTA file contents.

## Usage

In order to compare two initialized genomes one needs to issue the following command

```console
refgenie compare hg38 hg38_plus
```

The result is a set of informative flags that determine the level of compatibility of the genomes, for example:

```console
Flag: 2005
Binary: 0b11111010101

CONTENT_ALL_A_IN_B
LENGTHS_ALL_A_IN_B
NAMES_ALL_A_IN_B
CONTENT_A_ORDER
CONTENT_B_ORDER
CONTENT_ANY_SHARED
LENGTHS_ANY_SHARED
NAMES_ANY_SHARED
```

Based on the output above we can conclude that genome `hg38_plus` is a superset of `hg38`.
