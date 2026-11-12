# Which hg38 should I use?

The public server at <https://api.refgenie.org> holds about three dozen genomes related to hg38 (GRCh38). They share the same main chromosomes, but they differ in the extra pieces they carry, in how chromosomes are named, and in which pre-built assets they have. If you have no special need, pull `hg38`.

`hg38` is NCBI's "no-alt analysis set" (`GCA_000001405.15`). It has no ALT contigs, it has the most pre-built assets, and it is the version most short-read pipelines use.

Two facts make the choice easier:

1. **The main chromosomes match everywhere.** The GRCh38 genomes in the catalog all share the same lengths for chr1-22, X, Y and M. Only the names change (`chr1`, `1`, `NC_000001.11`).
2. **The differences that matter are small.** Genomes differ in which extra pieces they include (ALT contigs, decoys, HLA, EBV, unplaced pieces), whether a few regions are masked with N, and whether rare ambiguous bases (IUPAC codes such as R or M) are kept or turned into N. These change where reads land, so pick one genome and use it for every step.

## What to pull

| If you want to... | Pull | Why | Pre-built assets |
|---|---|---|---|
| Align short DNA reads and call variants (bwa, bowtie2) | `hg38` | No ALT contigs, so aligners that are not ALT-aware keep their mapping quality in duplicated regions. It includes chrEBV and the unplaced pieces. | FASTA and bwa, bowtie2, hisat2 and suffixerator indexes |
| The same, with decoys | `hg38_noalt_decoy` | `hg38` plus the hs38d1 decoy contigs, which absorb reads from DNA missing in the reference. Same coordinates as `hg38`. | The same as `hg38` |
| ALT-aware alignment (bwa-kit, GATK bundle) | `hs38DH` | `hg38` plus ALT contigs, decoys and HLA alleles. | None yet |
| RNA-seq with GENCODE annotation | `hg38` | GENCODE uses the same `chr1` ... `chrM` names and lengths, so its GTF lines up on the main chromosomes. GENCODE's own FASTA is `hg38_primary_gencode`. | `hg38` has a hisat2 index; `hg38_primary_gencode` has none |
| RNA-seq with Ensembl annotation (names `1`, `MT`) | `hg38_primary_ensembl` | Names match Ensembl GTFs. Same sequences as `hg38_primary_gencode`. | None |
| UCSC Genome Browser tracks | `hg38` | Every `hg38` name and length is also in UCSC's hg38, except chrEBV. UCSC's own file with ALT contigs is `hg38_initial`. | `hg38_initial` has FASTA and the same indexes as `hg38` |
| Only chromosome sizes, or one region | `hg38` | No need to download the genome. | `fasta_index`, or `refgenie getseq` |
| Only the 25 main chromosomes | `hg38_primary` | chr1-22, X, Y and M, and nothing else. | FASTA and bwa and bowtie2 indexes |

For example, to get the FASTA and a bwa index, or only the chromosome sizes, or one region:

```bash
refgenie pull hg38/fasta hg38/bwa_index
refgenie pull hg38/fasta_index
refgenie seek hg38/fasta_index.chrom_sizes
refgenie getseq -g hg38 -l chr1:10000-10100
```

Avoid for plain alignment any genome with ALT contigs or patches (`hg38_initial`, `hg38_p14`, `hg38_full_analysis`, `hg38_p14_genomic`, `hg38_p14_gencode`), unless your aligner is ALT-aware.

## How the hg38 genomes differ

"Masked" means chr5, 14, 19, 21, 22 and Y carry N in some regions. "IUPAC kept" means ambiguity codes are kept, not turned into N. Asset availability changes as builds finish, so check with `refgenie listr -g <alias>`.

| Alias | What it is | Sequences | ALT | Decoy | HLA | EBV | IUPAC kept | Masked |
|---|---|---|---|---|---|---|---|---|
| `hg38` | NCBI no-alt analysis set; same as `GRCh38-igenomes-ncbi` | 195 | - | - | - | yes | yes | yes |
| `hg38_analysis_set` | UCSC's copy of the same set; same as `hg38-igenomes-ucsc` | 195 | - | - | - | yes | no | yes |
| `hg38_noalt_decoy` | no-alt set plus hs38d1 decoys; same as `hg38_decoy_igenomes` | 2,580 | - | 2,385 | - | yes | yes | yes |
| `hg38_full_analysis` | NCBI full analysis set | 456 | 261 | - | - | yes | yes | yes |
| `hg38_full_analysis_plus_decoy` | full analysis set plus decoys | 2,841 | 261 | 2,385 | - | yes | yes | yes |
| `hs38DH` | 1000 Genomes / Broad GATK bundle | 3,366 | 261 | 2,385 | 525 | yes | yes | yes |
| `hg38_initial` | UCSC `hg38.fa`, 2013 release | 455 | 261 | - | - | - | no | - |
| `hg38_p14` | UCSC `hg38.p14.fa` | 711 | 261 + 256 patches | - | - | - | no | - |
| `hg38_primary` | UCSC main chromosomes only | 25 | - | - | - | - | no | - |
| `hg38_primary_gencode` | GENCODE / Ensembl 113 primary assembly; same as `hg38_primary_ensembl` | 194 | - | - | - | - | no | - |
| `hg38_p14_gencode` | GENCODE all / Ensembl 113 toplevel; same as `hg38_toplevel_ensembl` | 706 | 261 + 251 patches | - | - | - | no | - |
| `hg38_p14_genomic` | NCBI RefSeq, patch release 14 (`hg38_p0_genomic` is the first release) | 705 | 260-261 + patches | - | - | - | yes | - |
| `hg38_mm10` | `hg38_primary` plus `mm10_primary`, names prefixed `hg38_` and `mm10_` | 47 | - | - | - | - | no | - |

Every genome above except `hg38_primary` and `hg38_mm10` also has the 169 small unplaced and unlocalized pieces. Genomes marked "same as" hold the same sequences under different names or order. Soft-masked and unmasked files get the same digest, so they count as the same.

Names differ by source: UCSC and NCBI analysis sets use `chr1`, `chrM`, `chrUn_KI270302v1`; GENCODE uses `chr1`, `chrM`, `KI270302.1`; Ensembl uses `1`, `MT`, `KI270302.1`; RefSeq uses `NC_000001.11`. chrM is the same 16,569 bp sequence (rCRS) in every GRCh38 genome.

**`hg38` vs `hg38_primary`.** Both have the same 25 main chromosomes by name and length, but 18 of them differ in content: 12 because `hg38` keeps IUPAC codes, 4 because `hg38` is masked, and 2 for both reasons. `hg38` also has 170 extra pieces (169 unplaced and chrEBV).

To check how two genomes differ, use [`refgenie compare`](compare.md):

```bash
refgenie compare hg38 hg38_primary
```
