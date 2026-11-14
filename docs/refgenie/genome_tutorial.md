# Genome initialization and the RefgetStore

When you initialize a genome, refgenie loads its sequences into a local **RefgetStore**, a content-addressable store where each sequence is identified by its GA4GH digest. You can then retrieve any subsequence from the store, build a FASTA asset from it, and compare genomes.

## Prerequisites

- refgenie installed (`pip install "refgenie>=1.0.0a2"`; see [Install](install.md))
- refgenie initialized (`refgenie init`)
- the `fasta` recipe, from a synced [data channel](data_channels.md):

```bash
refgenie data-channel add refgenie https https://refgenie.github.io/refgenie-registry/index.yaml
refgenie data-channel sync refgenie --exists-ok
```

## Initialize a genome from a FASTA file

```bash
refgenie genome init \
  --fasta /path/to/hg38.fa.gz \
  --name hg38 \
  --description "Human genome build GRCh38" \
  --species "Homo sapiens"
```

Refgenie then:

1. Reads the FASTA file (`.fa`, `.fa.gz`, `.fasta`, or `.fasta.gz`).
2. Computes each sequence's GA4GH digest and stores the sequence in the RefgetStore at `~/.refgenie/genomes/.refget_store/`.
3. Computes a **sequence collection** digest for the whole genome.
4. Registers the genome under that digest, with the alias `hg38`.
5. Builds the genome's `fasta` asset. Use `--no-build` to skip this step.

Check that the genome is registered:

```bash
refgenie genome list
```

This shows each genome's digest, aliases, source, species, and description.

## Initialize a genome from a server

You can register a genome from a refgenie server or a RefgetStore. First, see what the server offers:

```bash
refgenie genome browse --server-url https://api.refgenie.org
```

Without `--server-url`, `genome browse` lists genomes on your subscribed servers. Then register one by its digest (`genome init --server` works for any genome the server lists). For example, to register hg38:

```bash
refgenie genome init \
  --server https://api.refgenie.org \
  --digest EiFob05aCWgVU_B_Ae0cypnQut3cxUP1 \
  --name hg38 \
  --no-build
```

Without `--no-build`, this downloads the genome's sequences to build its fasta asset. With `--no-build`, refgenie registers the genome without downloading them; `refgenie getseq` then reads only the bases you ask for. Use `--store <url>` instead of `--server` to register a genome straight from a RefgetStore. To register every genome on your subscribed servers, run `refgenie genome sync`.

## The fasta asset

The `fasta` asset is built from the sequences in the RefgetStore. It has three files:

| File | Seek key | Description |
|------|----------|-------------|
| `{digest}.fa` | `fasta` (default) | FASTA file with all sequences |
| `{digest}.fa.fai` | `fai` | FASTA index, computed by refgenie (no samtools needed) |
| `{digest}.chrom.sizes` | `chrom_sizes` | Sequence names and lengths |

For chromosome sizes alone, pull `fasta_index` (about 12 KB) and use `fasta_index.chrom_sizes`. It is the same file as `fasta.chrom_sizes`, without the 3 GB FASTA.

The files live under `~/.refgenie/genomes/data/{genome_digest}/fasta/{asset_digest}/`. `seek` returns paths in the alias tree, named for the alias you used:

```bash
refgenie seek hg38/fasta              # ~/.refgenie/genomes/alias/hg38/fasta/default/hg38.fa
refgenie seek hg38/fasta.fai          # .../hg38.fa.fai
refgenie seek hg38/fasta.chrom_sizes  # .../hg38.chrom.sizes
```

To rebuild the asset later, run `refgenie build hg38/fasta`.

## Retrieve sequences with getseq

`getseq` reads sequences straight from the RefgetStore, so you do not need a FASTA file on disk. It registers a genome from your servers on first use, so `genome init` is optional. With hg38 registered as above:

```bash
# Whole chromosome
refgenie getseq -g hg38 -l chrM

# Subsequence (0-based, half-open coordinates)
refgenie getseq -g hg38 -l chrM:0-60

# From a position to the end of the chromosome
refgenie getseq -g hg38 -l chrM:16500

# Write FASTA instead of plain bases
refgenie getseq -g hg38 -l chrM:0-60 --fasta
```

Refgenie resolves the alias `hg38` to its genome digest, looks up `chrM` in the sequence collection, and reads that sequence (or just the requested range) from the store. For a genome registered from a server, it reads only the bytes it needs. Add `--one-based` for samtools-style coordinates, for example `-l chrM:1-60`.

## Putting it all together

```bash
refgenie init
refgenie data-channel add refgenie https https://refgenie.github.io/refgenie-registry/index.yaml
refgenie data-channel sync refgenie --exists-ok

# Register a genome and build its fasta asset
refgenie genome init --fasta my_genome.fa.gz --name mygenome --description "My custom genome"
refgenie genome list

# Get paths to the fasta asset files
refgenie seek mygenome/fasta
refgenie seek mygenome/fasta.fai
refgenie seek mygenome/fasta.chrom_sizes

# Read sequence straight from the store
refgenie getseq -g mygenome -l chr1:0-500
```

## Key concepts

### The RefgetStore

The RefgetStore (`~/.refgenie/genomes/.refget_store/`) stores each sequence under its GA4GH digest. As a result:

- Identical sequences in different genomes are stored once.
- Any sequence can be checked against its digest.
- Subsequences can be read without loading a whole chromosome.

### Genome init and build

| | `genome init` | `build` |
|---|---|---|
| **Purpose** | Register a genome and load its sequences | Create asset files on disk |
| **Input** | FASTA file, server, or store | A registered genome (and any parent assets) |
| **Output** | Sequences in the store, a database record, and the `fasta` asset | Asset files, such as an aligner index |
| **Needed for** | `getseq`, `build`, `compare` | `seek` |

`getseq` works right after `genome init`.

## Compare genomes

Each genome is a GA4GH [sequence collection](https://ga4gh.github.io/seqcol-spec/): a set of sequences with names, lengths, and digests. Its identity is the sequence collection digest, so:

- The same sequences, with the same names and in the same order, give the same digest no matter who computes it or where.
- Adding or removing one sequence changes the digest.
- Renaming sequences (for example, `chr1` to `1`) changes the digest, but a comparison still shows that the sequences match.

`refgenie compare` reports how two genomes relate:

```bash
refgenie compare hg38 hg38_custom
```

The output lists which attributes (names, lengths, sequences, and their sorted forms) the two genomes share and how many elements each has. See [Compare genomes](compare.md).

## Next steps

Build more assets, stage them for serving, and push them to cloud storage with the [Building and serving tutorial](building_tutorial.md).
