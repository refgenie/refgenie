# Stage assets for serving

Staging gets a built asset ready to serve with `refgenie serve` or to push to cloud storage. An asset must be staged first.

## What staging does

Staging writes to the *stage folder*, in the form set by the asset class's [serving modes](serving_modes.md):

| Serving mode | What staging creates | Disk cost |
|---|---|---|
| `file` | A folder symlink to the built asset | None |
| `archive` | A `.tgz` archive | A full compressed copy |
| `none` | Nothing | None |

An asset class with both `file` and `archive` gets both.

## Stage, list, and unstage

```bash
refgenie stage add hg38/fasta
refgenie stage add hg38/fasta hg38/bowtie2_index mm10/fasta
refgenie stage add hg38/fasta --push-to my-s3-remote
refgenie stage list
refgenie stage remove hg38/fasta
```

`stage list` shows each staged asset's digest, registry path, mode, and archive size. `stage remove` deletes the symlink or archive; the built asset stays.

## Stage while building

```bash
refgenie build hg38/fasta --stage
refgenie build hg38/fasta --stage --push-to my-s3-remote
```

`--push-to` (on `build --stage` or `stage add`) marks the staged asset for upload to the named remotes, and `refgenie push` uploads it. `refgenie push` uploads only assets marked this way. See [Configure remote storage](remotes.md).

## Staging and serving

`refgenie serve` serves assets from the stage folder:

1. Build: `refgenie build hg38/fasta`
2. Stage: `refgenie stage add hg38/fasta`
3. Serve: `refgenie serve`

In `file` mode the server reads the built files through the symlink, so no data is copied. Cloud sync tools such as `aws s3 sync --follow-symlinks` can follow the symlinks too.

## The stage folder

Set the stage folder when you create the configuration:

```bash
refgenie init --genome-stage-folder /path/to/staged
```

If you leave it out, `refgenie init` uses `$REFGENIE_GENOME_STAGE_FOLDER`, or `~/.refgenie/archives` if that is not set. `refgenie config get` shows the current value.

For the whole build-to-serve workflow, see the [Building and serving tutorial](building_tutorial.md).
