# Upgrading from legacy refgenie

Refgenie 1.0 replaces the three legacy packages (`refgenie`, `refgenconf`, `refgenieserver`) with a single unified package. This tutorial walks you through migrating your existing data and configuration. For a reference of all API and CLI breaking changes, see the [migration reference](https://github.com/refgenie/refgenie1/blob/master/docs/migration-from-legacy.md).

!!! note "The upgrade tool is a separate package"
    The upgrade is performed by `refgenie-upgrade`, a standalone package. It is not a subcommand of `refgenie`.

## Before you start

You will need:

- Your legacy refgenie config file (the YAML file pointed to by `$REFGENIE`)
- Python 3.10 or higher
- A working [refgenie1 installation](install.md)

Back up your legacy config before starting:

```bash
cp $REFGENIE genome_config_backup.yaml
```

Install the upgrade tool:

```bash
pip install refgenie-upgrade
```

Because refgenie 1.0 is an alpha, install the prerelease first: `pip install --pre refgenie` (see the [install page](install.md)).

## Which path is yours?

**Path A: You built assets locally.** You ran `refgenie build`, have FASTA files on disk, or created custom assets. Follow the full upgrade below.

**Path B: You only pulled from a server.** You used `refgenie pull` to download pre-built assets from `refgenomes.databio.org` and have no local FASTAs or custom builds. Skip ahead to the [server consumer shortcut](#server-consumer-shortcut-path-b) at the bottom -- you just need a fresh `init` and `pull` from the new server.

## Preview what will be migrated (dry run)

Before making any changes, run a dry run to see what the tool will migrate:

```bash
refgenie-upgrade /path/to/genome_config.yaml --dry-run
```

The output lists the genomes, assets, and tags that will be migrated, along with any warnings about missing files or unsupported asset types. For example:

```
Summary: 2 genome(s) (5 asset(s)) would be migrated
```

Review the output before proceeding. If any genomes are skipped because their FASTA files are missing, see [Private genomes](#private-genomes-no-fasta-on-disk) in the troubleshooting section.

## Run the upgrade

```bash
refgenie-upgrade /path/to/genome_config.yaml
```

Or, if `$REFGENIE` is set:

```bash
refgenie-upgrade
```

This reads the legacy YAML config and:

1. Creates a new refgenie1 database (SQLite by default at `~/.refgenie/`)
2. Resolves genome digests to GA4GH format
3. Re-registers genome aliases (human-readable names like `hg38`)
4. Links existing asset files into the new directory layout
5. Creates alias symlinks for human-readable access

The original files are not moved or deleted. The upgrade creates new structures alongside them by default (using symlinks).

### Transfer modes

By default, `refgenie-upgrade` symlinks to your original files. You can change this:

```bash
refgenie-upgrade genome_config.yaml --symlink   # Default: symlink to originals
refgenie-upgrade genome_config.yaml --copy      # Copy files (doubles disk usage)
refgenie-upgrade genome_config.yaml --move      # Move files (frees legacy space)
```

### Custom target directory

To place the new refgenie1 instance in a different location:

```bash
refgenie-upgrade genome_config.yaml --target-dir /data/refgenie1
```

## Verify the upgrade

After the upgrade completes, verify that your assets are accessible:

```bash
refgenie1 list
```

You should see your genomes and assets listed. Check that you can retrieve paths:

```bash
refgenie1 seek hg38/fasta
refgenie1 seek hg38/bowtie2_index
```

The `seek` command now returns alias paths by default. Use `--abs` to get the digest-addressed content path. See [Retrieve paths to assets](seek.md) for details.

## Update your environment

Remove or rename the legacy `$REFGENIE` variable from your shell profile (`.bashrc`, `.zshrc`, etc.). If you leave it set, refgenie1 will print a warning each time it runs.

Here is how environment variables map from legacy to refgenie 1.0:

| Legacy | Refgenie 1.0 | Purpose |
|---|---|---|
| `$REFGENIE` | `$REFGENIE_DB_CONFIG_PATH` | Path to configuration |
| `genome_folder` (in YAML) | `$REFGENIE_GENOME_FOLDER` | Where assets are stored |
| `genome_archive_folder` (in YAML) | `$REFGENIE_GENOME_STAGE_FOLDER` | Where archives are staged |

Set the new variables only if you need non-default paths. See [Configuration](configuration.md) for the full list.

## Update pipeline configurations

### Seek-key changes

The `.dir` seek-key convention is removed. Asset classes now declare typed seek keys, and most assets have a single canonical seek key with the same name as the asset class. In practice, drop the `.dir` suffix:

```diff
- --genome-index { refgenie[sample.genome].bowtie2_index.dir }
+ --genome-index { refgenie[sample.genome].bowtie2_index }
```

Common changes:

| Legacy | Refgenie 1.0 |
|---|---|
| `hg38/bowtie2_index.dir` (returns directory; caller appends `/{genome}`) | `hg38/bowtie2_index` (returns prefix directly) |
| `hg38/bwa_index.dir` (caller appends `/{genome}.fa`) | `hg38/bwa_index` (returns prefix directly) |
| `hg38/fasta.dir` | `hg38/fasta` (returns the FASTA file path) |

For the full table of seek-key changes, see the [migration reference](https://github.com/refgenie/refgenie1/blob/master/docs/migration-from-legacy.md#seek-keys-the-dir-convention-is-gone).

## Server consumer shortcut (Path B)

If you only pulled pre-built assets and have no local FASTAs, start fresh instead of migrating:

```bash
pip install refgenie1
refgenie1 init
refgenie1 subscribe -s https://api.refgenie.org
refgenie1 pull hg38/fasta
```

No migration is needed. The server provides assets with GA4GH digests and the new directory layout.

## Troubleshooting

### Private genomes (no FASTA on disk)

Genomes are identified by GA4GH digests. If a genome is not in the built-in lookup table, the upgrade tool needs its FASTA on disk to compute the digest. If you built assets for a genome whose FASTA has been deleted, that genome will be skipped. To resolve this, re-download the FASTA and re-run the upgrade. You can also use the `--genome-folder` flag to point the tool at the correct location if your files have moved.

Note that legacy MD5 asset digests are not preserved; content digests are recomputed using SHA-256. Assets that used absolute paths (symlinked child assets) are also skipped.

### Partial failures

If the upgrade fails partway through (for example, from a full disk), re-run it. Already-migrated assets are skipped, so the process is safe to retry.

### Legacy $REFGENIE warning

If you see a warning about `$REFGENIE` being set, update your shell profile to remove or rename it. Refgenie 1.0 uses `$REFGENIE_DB_CONFIG_PATH` instead. See [Configuration](configuration.md) for details.

### What does not get migrated

The following are not carried over by the upgrade tool:

- Build logs and build recipes
- Archive configurations
- Server configurations
- Download statistics

These are internal bookkeeping data. If you need them, they remain in your legacy `genome_config.yaml` and genome folder.

## Cleanup

After you have verified that everything works:

1. Remove the legacy `genome_config.yaml` (you backed it up earlier).
2. Uninstall legacy packages: `pip uninstall refgenie refgenconf refgenieserver`
3. Remove `$REFGENIE` from your shell profile.

The legacy asset files under your genome folder remain in place. Refgenie 1.0 manages them through the new directory layout.
