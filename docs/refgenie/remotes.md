# Configure remote storage

Remotes are external storage endpoints (S3 buckets, HTTP servers) where you push staged assets for distribution. This page covers the client-side workflow for configuring remotes and pushing assets.

## Add a remote

Register a remote with a name, type, prefix URL, and optional push command:

```bash
refgenie remote add \
  --type s3 \
  --prefix my-bucket/refgenie \
  --name my-s3-remote \
  --push-command "aws s3 cp {local_path} s3://{prefix}/{relative_path}"
```

The `--type` is one of `http`, `https`, or `s3`. The `--prefix` is the base URL or path prefix where assets will be accessible after pushing. The `--name` identifies the remote in other commands (`push --remote`, `build --push-to`, `remote status`, `remote remove`); names must be unique and not all digits. The `--push-command` is an optional shell command template that refgenie executes when you run `refgenie push`.

### Push command placeholders

The push command template supports the following placeholders:

| Placeholder | Description |
|---|---|
| `{local_path}` | Absolute path to the local staged asset file |
| `{relative_path}` | Path of the asset relative to the stage folder |
| `{prefix}` | The prefix configured on the remote |
| `{genome_stage_folder}` | The full path to the genome stage folder (useful for folder sync strategies) |

**Per-asset push example** (uploads each asset individually):

```bash
--push-command "aws s3 cp {local_path} s3://{prefix}/{relative_path}"
```

**Folder sync push example** (syncs the entire stage folder):

```bash
--push-command "aws s3 sync {genome_stage_folder} s3://{prefix}/ --follow-symlinks"
```

## List remotes

View all configured remotes:

```bash
refgenie remote list
```

This displays each remote's ID, name, type, prefix, push command, and counts of pushed and unpushed assets.

## Remove a remote

Remove a remote by its name or ID:

```bash
refgenie remote remove my-s3-remote
```

## Check push status

See which assets have been pushed to each remote:

```bash
refgenie remote status
```

Filter by a specific remote:

```bash
refgenie remote status --remote my-remote-name
```

This shows per-remote counts of pushed and unpushed assets.

## Push workflow

The typical workflow for distributing assets to remote storage is:

1. **Build** the asset: `refgenie build hg38/fasta`
2. **Stage** the asset: `refgenie stage stage hg38/fasta`
3. **Push** assets: `refgenie push`

Push all unpushed assets to all configured remotes:

```bash
refgenie push
```

Push only assets for a specific genome:

```bash
refgenie push --genome hg38
```

Push to a specific remote only:

```bash
refgenie push --remote my-s3-remote
```

Preview what would be pushed without executing:

```bash
refgenie push --dry-run
```

You can also combine building, staging, and push intent in one step using the `build` command:

```bash
refgenie build hg38/fasta --stage --push-to my-s3-remote
refgenie push
```

When you push, refgenie executes the push command template configured on the remote and marks the asset as pushed.

### Push strategies

Refgenie supports two push strategies:

| Strategy | Description |
|---|---|
| `per_asset` (default) | Uploads each asset individually using `{local_path}` and `{relative_path}` |
| `folder_sync` | Syncs the entire `{genome_stage_folder}` to the remote |

Specify the strategy with `--strategy`:

```bash
refgenie push --strategy folder_sync
```

## Remote types

| Type | Use case |
|---|---|
| `http` | Assets served via HTTP (static hosting, CDN) |
| `https` | Assets served via HTTPS (S3 with HTTPS, CDN) |
| `s3` | Direct S3 integration |

## How remotes interact with the server

When a remote is configured on a refgenie server, the server redirects download requests to the remote URL instead of serving files directly. This allows the server to act as a metadata catalog while cloud storage handles the actual data transfer.
