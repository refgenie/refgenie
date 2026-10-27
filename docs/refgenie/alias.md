# Genome Aliases

## Introduction

Genome aliases provide human-readable names for reference genomes, mapping them to their underlying cryptographic digests. Instead of referring to genomes by their unique digest strings (e.g., `a5b35a4c5b27c7f3e8e2f4c3a1b2d9e7`), you can use memorable names like `hg38`, `mm10`, or `GRCh38`.

!!! success "Learning objectives"
    - What are genome aliases and why use them?
    - How to list existing aliases
    - How to set aliases for a genome
    - How to remove aliases

## Why use aliases?

Refgenie identifies genomes by their content-derived digest - a unique hash computed from the genome's sequences. While digests guarantee uniqueness and enable compatibility checking, they are difficult to remember and type. Aliases solve this problem by providing:

- **Human-readable names**: Use familiar names like `hg38` or `mm10` instead of digest strings
- **Multiple aliases per genome**: A single genome can have multiple aliases (e.g., `hg38`, `GRCh38`, `human`)
- **Consistent references**: Once set, aliases work throughout refgenie commands

!!! info "Digest vs Alias"
    A **digest** is the unique, content-derived identifier for a genome (computed from sequences). An **alias** is a human-friendly name that maps to a digest. You can use either in most refgenie commands.

## Listing aliases

To see all configured aliases, use the `alias get` command:

```console
refgenie alias get
```

This displays a table showing all aliases and their corresponding genome digests:

```console
┏━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Alias     ┃ Genome Digest                    ┃
┡━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ hg38      │ a5b35a4c5b27c7f3e8e2f4c3a1b2d9e7 │
│ GRCh38    │ a5b35a4c5b27c7f3e8e2f4c3a1b2d9e7 │
│ mm10      │ c3d45e6f7a8b9c0d1e2f3a4b5c6d7e8f │
└───────────┴──────────────────────────────────┘
```

### Filtering by alias

To look up the digest for specific aliases:

```console
refgenie alias get -a hg38 GRCh38
```

### Filtering by digest

To find all aliases for a specific genome digest:

```console
refgenie alias get -g a5b35a4c5b27c7f3e8e2f4c3a1b2d9e7
```

!!! note "Mutually exclusive options"
    The `-a/--aliases` and `-g/--genome-digests` options are mutually exclusive. Use one or the other to filter results.

## Setting aliases

To assign an alias to a genome, use the `alias set` command with the `-a` (alias) and `-d` (digest) options:

```console
refgenie alias set -a hg38 -d a5b35a4c5b27c7f3e8e2f4c3a1b2d9e7
```

### Setting multiple aliases at once

You can set multiple aliases for the same genome in a single command:

```console
refgenie alias set -a hg38 GRCh38 human -d a5b35a4c5b27c7f3e8e2f4c3a1b2d9e7
```

### Replacing existing aliases

The `--reset` flag removes all existing aliases for a genome before setting new ones:

```console
refgenie alias set -a hg38_v2 -d a5b35a4c5b27c7f3e8e2f4c3a1b2d9e7 --reset
```

This is useful when you want to rename a genome or consolidate aliases. Without `--reset`, new aliases are added to any existing ones.

!!! warning "Genome must exist"
    You can only set aliases for genomes that already exist in your refgenie database. If you try to set an alias for a non-existent digest, you will see an error:

    ```console
    Genome with digest abc123... does not exist. You must initialize it first by building/pulling an asset for that genome.
    ```

## Removing aliases

To remove aliases, use the `alias remove` command:

```console
refgenie alias remove -a hg38
```

### Removing multiple aliases

You can remove multiple aliases at once:

```console
refgenie alias remove -a hg38 GRCh38
```

!!! tip "Aliases vs Genome"
    Removing an alias does not delete the genome or its assets. It only removes the name mapping. The genome remains accessible by its digest or any remaining aliases.

## Using aliases in commands

Once set, aliases can be used anywhere a genome identifier is expected:

```console
# Pull an asset using an alias
refgenie pull hg38/fasta

# List assets for an aliased genome
refgenie list -g hg38

# Get sequence using an alias
refgenie getseq -g hg38 -l chr1:1-100

# Compare genomes using aliases
refgenie compare hg38 mm10
```

## Command reference

### alias get

List aliases with optional filtering.

```console
refgenie alias get [-a ALIASES...] [-g DIGESTS...]
```

| Option | Description |
|--------|-------------|
| `-a, --aliases` | Filter by specific alias names |
| `-g, --genome-digests` | Filter by specific genome digests |

### alias set

Set one or more aliases for a genome digest.

```console
refgenie alias set -a ALIASES... -d DIGEST [--reset]
```

| Option | Description |
|--------|-------------|
| `-a, --aliases` | One or more alias names to set (required) |
| `-d, --digest` | Genome digest to associate with the aliases |
| `-r, --reset` | Remove all existing aliases before setting new ones |
| `-f, --force` | Force the action if genome does not exist |

### alias remove

Remove one or more aliases.

```console
refgenie alias remove -a ALIASES...
```

| Option | Description |
|--------|-------------|
| `-a, --aliases` | One or more alias names to remove (required) |
