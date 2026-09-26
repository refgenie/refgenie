# Stage assets for serving

Staging prepares built assets for serving to clients or pushing to cloud storage. Before an asset can be served or pushed, it must be staged.

## What staging does

Staging creates entries in the `genome_stage_folder` (default: `~/.refgenie/archives/`) based on the asset's serving modes. The serving mode is defined by the asset class -- see [Serving modes](serving_modes.md) for background.

| Serving mode | What staging creates | Disk cost |
|---|---|---|
| `file` | A directory symlink from the stage folder to the asset's build directory | None (symlink only) |
| `archive` | A `.tgz` tarball in the stage folder | Full compressed copy |
| `none` | Nothing -- staging is skipped | None |

If an asset class specifies both `file` and `archive` modes, staging creates both a symlink and a tarball.

## Stage an asset

After building an asset, stage it:

```bash
refgenie stage stage hg38/fasta
```

You can stage multiple assets at once:

```bash
refgenie stage stage hg38/fasta hg38/bowtie2_index mm10/fasta
```

## List staged assets

View all staged assets with their mode and size:

```bash
refgenie stage staged
```

This displays a table with the asset digest, registry path, mode (file or archive), and size (for archive mode).

## Unstage an asset

Remove staging records and clean up staged files:

```bash
refgenie stage unstage hg38/fasta
```

## Build, stage, and push in one step

The `build` command supports `--stage` and `--push-to` flags to combine building, staging, and push intent creation:

```bash
refgenie build hg38/fasta --stage --push-to my-s3-remote
```

The `--stage` flag stages the asset after building. The `--push-to` flag creates push intent records for the specified remote(s) after staging. See [Configure remote storage](remotes.md) for setting up remotes and pushing.

## How staging relates to serving

The refgenie server (`refgenie serve`) serves assets from the stage folder. The workflow is:

1. **Build** the asset: `refgenie build hg38/fasta`
2. **Stage** the asset: `refgenie stage stage hg38/fasta`
3. **Serve** the asset: `refgenie serve`

File-mode staging is especially efficient for local servers because the symlink means no data is duplicated -- the server reads directly from the build directory through the symlink. Cloud sync tools like `aws s3 sync` also follow symlinks by default.

For the complete build-to-serve workflow, see the [Building and serving tutorial](building_tutorial.md).

## Stage folder configuration

The stage folder location is controlled by the `REFGENIE_GENOME_STAGE_FOLDER` environment variable:

```bash
export REFGENIE_GENOME_STAGE_FOLDER=~/.refgenie/archives
```

The default is `$REFGENIE_HOME_PATH/archives`.
