# Which mm39 should I use?

The public server at <https://api.refgenie.org> holds about a dozen genomes for mouse mm39 (GRCm39). Only one of them has pre-built assets: `mm39`. If you have no special need, pull `mm39`.

`mm39` is the GRCm39 primary assembly from Ensembl release 113 (`GCA_000001635.9`). It has 61 sequences: chromosomes `1` to `19`, `X`, `Y` and `MT`, plus 39 unplaced and unlocalized pieces.

**Watch the names.** `mm39` uses Ensembl names (`1`, `MT`). Refgenie's `hg38` and `mm10` use UCSC names (`chr1`, `chrM`). A BED or GTF file with `chr1` names will not line up with reads aligned to `mm39` until you rename one of them.

## What to pull

| If you want to... | Pull | Notes |
|---|---|---|
| Align reads (bwa, bowtie2, hisat2) | `mm39` | Has FASTA and bwa, bowtie2, hisat2 and suffixerator indexes. |
| RNA-seq with Ensembl annotation | `mm39` | Names match Ensembl GTFs. |
| Work with `chr1`-style names (UCSC, GENCODE) | `mm39`, then rename | No chr-named mm39 has pre-built assets yet. See below. |
| Only chromosome sizes | `mm39` | Pull `fasta_index`, about 3 KB. |
| The older mm10 (GRCm38) | `mm10` | NCBI's no-alt analysis set, with `chr1` names and the same indexes as `mm39`. |

```bash
refgenie pull mm39/bowtie2_index
refgenie pull mm39/fasta_index
refgenie seek mm39/fasta_index.chrom_sizes
```

## If your files use chr1-style names

These genomes hold the same 61 sequences as `mm39`, only with other names and in another order, so the coordinates line up once the names match:

- `mm39_ucsc`: UCSC names, such as `chr1`, `chrM` and `chr1_GL456210v1_random`.
- `mm39_gencode`: GENCODE names, such as `chr1`, `chrM` and `GL456210.1`.
- `mm39_full_analysis`: NCBI's analysis set, with UCSC-style names.

None of them has pre-built assets yet. You can:

1. Use `mm39` and rename the chromosomes in your annotation (`chr1` to `1`, `chrM` to `MT`). `refgenie compare mm39 mm39_ucsc` prints every rename.
2. Build the assets yourself from your FASTA file. See [Build assets](build.md).
3. Ask for them with the "Request assets" link on the genome's page at [refgenie.org](https://refgenie.org).

Avoid the masked copies (`mm39_ucsc_masked`, `mm39_primary_rm_ensembl`, `mm39_toplevel_rm_ensembl`) for alignment: repeats in them are replaced with N.

For human genomes, see [Which hg38 should I use?](choosing_hg38.md). For yeast, see [Which yeast genome should I use?](choosing_yeast.md).
