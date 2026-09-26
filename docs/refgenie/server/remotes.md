# Understanding Remotes in Refgenieserver

Refgenieserver can be configured to work with "remotes," which are external storage locations for asset archives. This allows the server to provide access to assets that may not be stored directly on the server's local filesystem, offering flexibility in how and where large asset files are hosted.

Examples of remote types include: `http`, `https`, `s3`

## Serving Files: Remote Redirect or Local Serve via `/archives/{digest}/download`

The `/archives/{digest}/download` route is responsible for serving asset archive files. The server determines how to serve the file based on the following logic:

### 1. Attempt Remote Redirect (HTTP/HTTPS)

- When a request comes in for an asset archive (identified by its `digest`):
  - Refgenieserver first retrieves the archive's metadata, including its local path (`archive.archive_path`).
  - It then checks if an HTTP or HTTPS remote is configured and if the archive's local path can be resolved to a URL on that remote.
    - This involves querying the database for `Remote` configurations of type `http` or `https` and the server's `genome_archive_folder` (from the `Configuration` table). HTTPS remotes are prioritized.
    - If a suitable remote is found, and the archive's path is relative to the `genome_archive_folder`, a full remote URL is constructed (e.g., `https://remote.example.com/path/to/archive.tgz`).
  - If a remote URL is successfully constructed, the server issues an HTTP redirect to this URL. The client then fetches the archive directly from the remote server.

### 2. Local File Serve (Fallback)

- If no suitable HTTP/HTTPS remote is configured, or if the archive's path cannot be resolved to a remote URL (e.g., it's not relative to the configured `genome_archive_folder`), Refgenieserver falls back to serving the file directly from its local filesystem.
- The server uses the `archive.archive_path` to locate the file on its disk.
- If the file exists locally, it is streamed to the client as a direct download.
- If the archive metadata or the local file itself is not found, an appropriate HTTP error (e.g., 404 Not Found) is returned.

> This mechanism ensures that production Refgenieserver instances can serve files from blob storage, while also enabling local development and testing by allowing local file serving when no remotes are configured.

## Managing Remotes

Remotes are managed using the `refgenie remote` command-line interface.

### Constraint: One Remote Per Type

A crucial rule for managing remotes is that **you can only have one remote of each type configured at any given time**. For example, you can have one S3 remote and one HTTP remote, but you cannot have two S3 remotes. If you try to add a remote of a type that already exists, the new configuration will typically overwrite the old one for that type.

### Adding a Remote

To add a new remote, you use the `refgenie remote add` command. You need to specify a name, the type of the remote, and a prefix (the base path or URL).

**Command:**

```bash
refgenie remote add --type <remote_type> --prefix <prefix_value> --name <remote_name>
```

**Arguments:**

- `--type`: (Required) The type of the remote. Examples: `http`, `s3`. This is chosen from the available `RemoteType` values in Refgenie.
- `--prefix`: (Required) The base identifier for the remote.
  - For `http`, this would be the base URL (e.g., `https://example.com/assets/`).
  - For `s3`, this would be the S3 bucket name (e.g., `s3://my-refgenie-bucket/`).
- `--name`: (Required) The name you use to refer to the remote in other commands. Names must be unique and not all digits.

**Example (adding an S3 remote):**

```bash
refgenie remote add --type s3 --prefix "s3://my-assets-bucket/refgenie_archives/" --name primary-s3
```

### Listing Configured Remotes

To see the currently configured remotes, use the `refgenie remote list` command.

**Command:**

```bash
refgenie remote list
```

This will output a table listing all configured remotes, their names, types, and prefixes. 

### Removing a Remote

To remove an existing remote, you use the `refgenie remote remove` command, giving the name or ID of the remote you wish to delete.

**Command:**

```bash
refgenie remote remove <remote_name_or_id>
```

**Arguments:**

- `<remote_name_or_id>`: (Required) The name or numeric ID of the remote to remove.

**Example (removing an S3 remote):**

```bash
refgenie remote remove primary-s3
```

This removes that remote from your Refgenie setup. If no remote has that name or ID, the command exits with a "not found" error.
