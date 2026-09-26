# Genome initialization and the RefgetStore

This tutorial walks through the new refgenie workflow for initializing genomes, building FASTA assets, and retrieving sequences. The key concept: when you initialize a genome, refgenie loads all sequences into a local **RefgetStore** -- a content-addressable store where every sequence is identified by its GA4GH digest. This means you can retrieve any subsequence without needing FASTA files on disk.

## Prerequisites

- refgenie installed (`pip install refgenie1-*.whl`)
- refgenie initialized (`refgenie1 init`)

## Initialize a genome from a local FASTA file

The `genome init` command registers a genome and loads its sequences into the RefgetStore:

```bash
refgenie1 genome init \
  --fasta /path/to/hg38.fa.gz \
  --name hg38 \
  --description "Human genome build GRCh38" \
  --species "Homo sapiens"
```

What happens under the hood:

1. Refgenie reads the FASTA file (supports `.fa`, `.fa.gz`, `.fasta`, `.fasta.gz`)
2. Each sequence is digested (SHA-512/24 truncated) and stored in the RefgetStore at `~/.refgenie/genomes/.refget_store/`
3. A **sequence collection** digest is computed from the set of all sequences
4. The genome is registered in the database with that digest, along with the alias `hg38`

You can verify the genome was added:

```bash
refgenie1 genome list
```

This shows the genome digest, aliases, source (local), and description.

## Initialize a genome from a remote seqcolapi server

You can also register a genome from a remote Sequence Collections API server without downloading any sequence data:

```bash
refgenie1 genome init \
  --remote-url https://seqcolapi.databio.org \
  --digest abc123def456... \
  --name hg38_remote \
  --description "Human GRCh38 (remote)"
```

This verifies the collection exists on the remote server and stores only the metadata locally. No sequences are downloaded -- they're fetched on demand when you use `getseq`.

!!! tip "Browsing remote genomes"
    Use `genome browse` to see what's available on a remote server:
    ```bash
    refgenie1 genome browse https://seqcolapi.databio.org
    ```

## Build a FASTA asset from the store

Once a genome is initialized, you can build a FASTA asset. This **exports** sequences from the RefgetStore into standard FASTA files:

```bash
refgenie1 build hg38/fasta
```

This creates three files:

| File | Description |
|------|-------------|
| `{digest}.fa` | FASTA file with all sequences |
| `{digest}.fa.fai` | FASTA index (computed natively, no samtools needed) |
| `{digest}.chrom.sizes` | Chromosome names and lengths |

The files are stored under `~/.refgenie/genomes/data/{genome_digest}/fasta/default/`.

!!! note "No samtools required"
    The FASTA index is computed by refgenie itself, so you don't need samtools installed for the `fasta` recipe.

## Seek asset paths

Use `seek` to get the path to any file in the built asset:

```bash
# Default seek key (the FASTA file)
refgenie1 seek hg38/fasta

# Specific seek keys
refgenie1 seek hg38/fasta.fasta      # → /path/to/{digest}.fa
refgenie1 seek hg38/fasta.fai        # → /path/to/{digest}.fa.fai
refgenie1 seek hg38/fasta.chrom_sizes # → /path/to/{digest}.chrom.sizes
```

Use these in your pipelines instead of hardcoded paths:

```bash
bwa mem $(refgenie1 seek hg38/fasta) reads.fq > aligned.sam
```

## Retrieve sequences with getseq

The `getseq` command retrieves sequences directly from the RefgetStore -- no FASTA file needed on disk. This is the refget-style sequence retrieval:

```bash
# Whole chromosome
refgenie1 getseq -g hg38 -l chr1

# Subsequence (0-based, half-open coordinates)
refgenie1 getseq -g hg38 -l chr1:0-1000

# From a position to end of chromosome
refgenie1 getseq -g hg38 -l chr1:50000
```

### How it works

1. The alias `hg38` is resolved to its genome digest
2. The sequence name (`chr1`) is looked up in the sequence collection
3. The individual sequence digest is used to retrieve data from the RefgetStore
4. For subsequences, only the requested range is extracted

### Remote sequence retrieval

If the genome was initialized from a remote server, `getseq` fetches sequences on demand and caches them locally:

```bash
# First call: fetches chr1 from the remote server and caches it
refgenie1 getseq -g hg38_remote -l chr1:0-100

# Subsequent calls: served from the local cache
refgenie1 getseq -g hg38_remote -l chr1:500-600
```

## Putting it all together

Here's a complete workflow from FASTA file to sequence retrieval:

```bash
# 1. Initialize refgenie
refgenie1 init

# 2. Initialize a genome from a local FASTA
refgenie1 genome init \
  --fasta my_genome.fa.gz \
  --name mygenome \
  --description "My custom genome"

# 3. Verify it's registered
refgenie1 genome list

# 4. Build the FASTA asset (exports from store to files)
refgenie1 build mygenome/fasta

# 5. Seek paths to the built files
refgenie1 seek mygenome/fasta           # FASTA file
refgenie1 seek mygenome/fasta.fai       # FASTA index
refgenie1 seek mygenome/fasta.chrom_sizes  # Chrom sizes

# 6. Retrieve sequences directly from the store
refgenie1 getseq -g mygenome -l chr1:0-500
```

## Key concepts

### The RefgetStore

The RefgetStore (located at `~/.refgenie/genomes/.refget_store/`) is a content-addressable database of sequences. Every sequence is identified by its GA4GH SHA-512/24 digest. This means:

- **Deduplication**: Identical sequences across genomes are stored only once
- **Verification**: Sequences can be verified against their digest at any time
- **Direct access**: Subsequences can be retrieved without loading entire chromosomes

### Genome init vs. build

These serve different purposes:

| | `genome init` | `build` |
|---|---|---|
| **Purpose** | Register a genome, load sequences into RefgetStore | Create asset files on disk |
| **Input** | FASTA file or remote URL | Sequences already in RefgetStore |
| **Output** | Sequences in store + database record | Files on disk (`.fa`, `.fai`, `.chrom.sizes`) |
| **Required for** | `getseq`, `build`, `compare` | `seek` (to get file paths for tools) |

You must `genome init` before you can `build`. But you can use `getseq` immediately after `genome init` -- no `build` needed.

## Sequence collections and genome comparison

The RefgetStore uses the GA4GH **sequence collections** (seqcol) standard to represent genomes. Each genome is a sequence collection -- a set of sequences with names, lengths, and GA4GH digests. This enables powerful comparison features.

### Compare two genomes

Use `refgenie compare` to see how two genomes relate to each other:

```bash
refgenie compare hg38 hg38_custom
```

This performs a sequence collection comparison and reports:
- Which sequences are shared between the two genomes
- Which sequences are unique to each genome
- Whether the genomes differ only in sequence names, ordering, or actual sequence content

### How sequence collection digests work

Each genome's identity is a **sequence collection digest** -- computed from the sorted, digested set of all its sequences. This means:

- Two genomes built from the same FASTA file will always have the same digest, regardless of who built them or where
- Adding or removing even one sequence changes the digest
- Renaming sequences (e.g., `chr1` vs `1`) produces a different digest, but `compare` can still detect that the underlying sequences are identical

### Integration with remote seqcol servers

The RefgetStore is compatible with the [GA4GH seqcol specification](https://ga4gh.github.io/seqcol-spec/). You can initialize genomes from remote seqcol servers and compare local genomes against remote collections:

```bash
# Initialize from a remote seqcol server
refgenie genome init \
  --store https://seqcolapi.databio.org \
  --digest abc123... \
  --name hg38_remote

# Browse available genomes on a remote server
refgenie genome browse https://seqcolapi.databio.org
```

This interoperability means refgenie can participate in the broader GA4GH ecosystem for genome identification and comparison.

## Next steps

Now that you have a genome initialized, you can build assets, stage them for serving, and push them to cloud storage. The [Building and Serving Tutorial](building_tutorial.md) walks through the complete lifecycle step by step.
