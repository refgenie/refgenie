# Upgrade from legacy refgenie

This guide moves a refgenie 0.x setup (config versions 0.3 and 0.4) to refgenie 1.0. The `refgenie-upgrade` tool reads your old config, registers your genomes and assets in a new 1.0 database, and links your existing asset files, so you do not have to rebuild anything. For the changes that affect code and pipelines, see the [migration reference](migration-from-legacy.md).

## Before you start

You need:

- Your 0.x config file (the YAML file that `$REFGENIE` points to)
- Python 3.11 or higher

Back up the config:

```bash
cp $REFGENIE genome_config_backup.yaml
```

Install the upgrade tool from GitHub. It also installs refgenie 1.0 if you do not have it.

```bash
pip install git+https://github.com/refgenie/refgenie-upgrade.git
```

## Which path is yours?

- **You built assets locally** (you ran `refgenie build`, have FASTA files on disk, or added custom assets): follow the steps below.
- **You only pulled assets from a server**: skip to [Start fresh](#start-fresh-if-you-only-pulled-assets).

## Step 1: Preview the upgrade

A dry run shows what the tool will migrate, without changing anything:

```bash
refgenie-upgrade /path/to/genome_config.yaml --dry-run
```

It lists the genomes, assets, and tags it will migrate, warns about missing files, and ends with a summary such as:

```
Summary: 2 genome(s) (5 asset(s)) would be migrated
```

If a genome is skipped because its FASTA file is missing, see [Troubleshooting](#troubleshooting).

## Step 2: Run the upgrade

```bash
refgenie-upgrade /path/to/genome_config.yaml
```

If `$REFGENIE` is set, you can leave out the path. The tool:

1. Creates a new refgenie 1.0 instance (database and genome folder).
2. Computes GA4GH digests for your genomes.
3. Registers your genome aliases, like `hg38`.
4. Links your existing asset files into the new folder layout.

By default the new instance goes in a `refgenie1_data` folder next to your old genome folder. Use `--target-dir` to put it somewhere else:

```bash
refgenie-upgrade genome_config.yaml --target-dir /data/refgenie
```

By default the tool symlinks to your original files and leaves them in place. You can copy or move them instead:

```bash
refgenie-upgrade genome_config.yaml --symlink   # default: symlink to originals
refgenie-upgrade genome_config.yaml --copy      # copy files (doubles disk use)
refgenie-upgrade genome_config.yaml --move      # move files (frees the old space)
```

## Step 3: Point refgenie at the new instance

When it finishes, the tool prints the line to add to your shell profile (`.bashrc`, `.zshrc`, etc.), for example:

```bash
export REFGENIE_HOME_PATH=/data/refgenie1_data
```

Also remove `$REFGENIE` from your shell profile. Refgenie 1.0 does not read it and prints a warning while it is set.

| 0.x setting | 1.0 setting | Purpose |
|---|---|---|
| `$REFGENIE` | `$REFGENIE_HOME_PATH` (or `$REFGENIE_DB_CONFIG_PATH`) | Where refgenie finds its configuration |
| `genome_folder` (in YAML) | `$REFGENIE_GENOME_FOLDER` | Where assets are stored |
| `genome_archive_folder` (in YAML) | `$REFGENIE_GENOME_STAGE_FOLDER` | Where assets are staged for serving |

See [Configure](configuration.md) for all settings.

## Step 4: Check the result

```bash
refgenie list
refgenie seek hg38/fasta
refgenie seek hg38/bowtie2_index
```

`seek` returns paths in the human-readable alias tree. Add `--abs` for the digest-addressed path. See [Retrieve paths to assets](seek.md).

## Step 5: Update your pipelines

Seek keys ending in `.dir` no longer exist. In most cases, drop the suffix; aligner indexes now return the index prefix directly:

```diff
- --genome-index { refgenie[sample.genome].bowtie2_index.dir }
+ --genome-index { refgenie[sample.genome].bowtie2_index }
```

The [migration reference](migration-from-legacy.md) lists this and the other changes to seek keys, the Python API, looper integration, and plugins.

## Start fresh if you only pulled assets

If you only pulled pre-built assets, set up a new instance and pull again from the public server:

```bash
pip install "refgenie>=1.0.0a1"
refgenie init
refgenie pull hg38/fasta
```

`refgenie init` subscribes to `https://api.refgenie.org`, which serves assets in the 1.0 format.

## Troubleshooting

**A genome was skipped because its FASTA is missing.** The tool knows the digests of common genomes. For other genomes it needs the FASTA file on disk to compute the digest. Put the FASTA back and run the upgrade again. If your files have moved, point the tool at them with `--genome-folder`.

**The upgrade stopped partway**, for example because the disk filled up. Run it again. It skips assets that are already migrated.

**Some things are not migrated.** The tool skips assets stored at absolute paths (symlinked child assets). It recomputes asset digests with SHA-256 instead of keeping the old MD5 digests. Build logs, build recipes, archive and server settings, and download statistics stay in your old config file and genome folder.

## Clean up

After everything works:

1. Delete the old `genome_config.yaml` (you have a backup).
2. Uninstall the old packages: `pip uninstall refgenconf refgenieserver`.
3. Make sure `$REFGENIE` is gone from your shell profile.

If you used the default `--symlink` mode, keep your old genome folder: the new instance links to the files in it.
