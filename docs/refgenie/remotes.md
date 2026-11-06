# Configure remote storage

A *remote* is cloud storage or a web host (such as an S3 bucket) where you upload staged assets so others can download them. You add a remote, mark staged assets for it, and run `refgenie push`.

## Add a remote

```bash
refgenie remote add \
  --type s3 \
  --prefix my-bucket/refgenie \
  --name my-s3-remote \
  --push-command "aws s3 cp {local_path} s3://{prefix}/{relative_path}"
```

- `--type` is `http`, `https`, or `s3`.
- `--prefix` is the base address where the files will be found after upload.
- `--name` is what other commands use to refer to the remote (`build --push-to`, `stage add --push-to`, `push --remote`, `remote status`, `remote remove`). It must be unique and not all digits.
- `--push-command` is the shell command refgenie runs to upload. A remote without one is skipped by `refgenie push`.

The push command can use these placeholders:

| Placeholder | Value |
|---|---|
| `{local_path}` | Full path to the staged file or folder |
| `{relative_path}` | Its path inside the stage folder |
| `{prefix}` | The remote's prefix |
| `{genome_stage_folder}` | The stage folder (for the `folder_sync` strategy) |

## Push

`refgenie push` uploads assets that were staged with `--push-to`:

```bash
refgenie build hg38/fasta --stage --push-to my-s3-remote
refgenie stage add mm10/fasta --push-to my-s3-remote  # an asset already built
refgenie push
```

Narrow it down, or preview it:

```bash
refgenie push --genome hg38          # one genome
refgenie push --remote my-s3-remote  # one remote
refgenie push --dry-run              # show what would run
```

Refgenie runs the remote's push command for each asset and records it as pushed. For `s3` remotes it checks that the upload landed, and skips archives that are already in the bucket.

### Sync the whole stage folder

The default strategy, `per_asset`, runs the push command once per asset. The `folder_sync` strategy runs it once per remote, so the command should copy the whole stage folder:

```bash
refgenie remote add \
  --type s3 \
  --prefix my-bucket/refgenie \
  --name my-s3-sync \
  --push-command "aws s3 sync {genome_stage_folder} s3://{prefix}/ --follow-symlinks"
refgenie push --strategy folder_sync
```

Refgenie then marks each waiting asset as pushed if its staged file exists.

## List, check, and remove remotes

```bash
refgenie remote list                         # each remote and its pushed/unpushed counts
refgenie remote status                       # pushed/unpushed counts and unpushed digests
refgenie remote status --remote my-s3-remote
refgenie remote remove my-s3-remote          # by name or id
```

## Remotes and the server

When a refgenie server knows an asset was pushed to a remote, it sends download requests to the remote's address instead of serving the file itself. The server then acts as a catalog while the cloud storage moves the data.
