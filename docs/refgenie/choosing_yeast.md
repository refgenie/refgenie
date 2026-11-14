# Which yeast genome should I use?

For budding yeast (*Saccharomyces cerevisiae*, strain S288C), two genomes on the public server have pre-built assets: `R64-1-1` and `sacCer3`. Both are the same assembly, R64-1-1, with the same 17 sequences in the same order. Only the names differ. Pick the one whose names match your other files.

| Pull | Names | Pre-built assets |
|---|---|---|
| `R64-1-1` (also called `yeast_s288c`) | `I` to `XVI`, and `Mito`. From Ensembl release 81. | FASTA and bwa, bowtie2, hisat2 and suffixerator indexes |
| `sacCer3` | `chrI` to `chrXVI`, and `chrM`. UCSC's names, as packaged by Illumina iGenomes. | FASTA and bwa, bowtie2 and hisat2 indexes |

```bash
refgenie pull R64-1-1/bowtie2_index
refgenie compare sacCer3 R64-1-1   # prints every rename: chrI -> I, ..., chrM -> Mito
```

## sacCer3 from UCSC is not the same file

Refgenie's `sacCer3` is the iGenomes copy, which lists `chrM` sixth, after `chrIX`. The `sacCer3.fa.gz` you download from UCSC has the same 17 sequences with the same names, but lists `chrM` last. Sequence order is part of a genome's digest, so UCSC's file has a different digest (`nHYYtlwZOhq203qjGUFNUaJlgiO2ISMs`), and refgenie does not list it. If you mix the two, your BAM files and your reference list the chromosomes in different orders, which some tools, such as GATK, reject.

Ensembl's release 81 file, `Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz`, has exactly the same digest as refgenie's `R64-1-1`.

To check a FASTA file you already have, compute its digest with `refget`, which comes with refgenie, and look it up:

```bash
refget fasta digest sacCer3.fa.gz
curl -s https://api.refgenie.org/v4/genomes/<digest>
```

A `404` means the server does not have that exact genome.

Older yeast assemblies, such as `sacCer2-igenomes-ucsc`, are listed too, but have no pre-built assets.

For human and mouse genomes, see [Which hg38 should I use?](choosing_hg38.md) and [Which mm39 should I use?](choosing_mm39.md).
