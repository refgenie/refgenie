# Remotes on a server

A remote is external storage, such as an S3 bucket or a web server, that holds copies of your staged assets. When a server has pushed an asset to a remote, it can send downloads there instead of serving the file itself. To add remotes and push assets, see [Configure remote storage](../remotes.md).

## How the server answers a download

For `/v4/archives/{asset_digest}/download` and `/v4/assets/{asset_digest}/files/{file_path}`:

1. **Redirect.** If the asset was pushed to an `http` or `https` remote, the server redirects to it. The URL is the remote's prefix plus the file's path inside the stage folder, so a staged archive at `<stage folder>/<genome digest>/fasta/<asset digest>.tgz` with prefix `https://example.org/refgenie` is sent to `https://example.org/refgenie/<genome digest>/fasta/<asset digest>.tgz`. When several remotes qualify, `https` wins over `http`, then the newest one.
2. **Serve locally.** Otherwise the server streams the file from its stage folder, or returns 404 if the file is missing.

Only `http` and `https` remotes give a download URL; an `s3` remote is for pushing. If you push to S3 and serve from a separate server, `refgenie catalog-export --https-prefix <public bucket URL>` writes a catalog whose download links point at the bucket's public URL. See the [CLI reference](../refgenie-cli.md#catalog-export).

## Manage remotes

```bash
refgenie remote add --type https --prefix https://example.org/refgenie --name public
refgenie remote list
refgenie remote remove public
```

You can have several remotes of each type. Names must be unique and not all digits. See [Configure remote storage](../remotes.md) for push commands and push status.
