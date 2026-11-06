# Serving modes

An asset class's seek keys say what an asset contains. Its `serving_modes` field says how a server delivers it. A FASTA file is often fetched on its own, while a bowtie2 index is six files that are only useful together, so the two are best delivered in different ways.

To see serving modes in practice, including staging and serving file-mode and archive-mode assets, see the [Building and serving tutorial](building_tutorial.md).

## The three modes

| Mode | Server stores | On request |
|------|--------------|------------|
| `file` | The asset's individual files | Serves each file at its own URL |
| `archive` | A `.tgz` tarball | Serves the tarball for bulk download |
| `none` | Nothing but metadata | Returns metadata, digests and build commands |

An asset class lists one or more modes:

```yaml
serving_modes:
  - file
```

When an asset class does not set `serving_modes`, it defaults to `[file]`.

- **File mode** lets clients fetch one file with a plain HTTP request, without downloading the whole asset. It suits web apps, cloud pipelines and other HTTP clients.
- **Archive mode** sends the whole asset as one compressed tarball, which the client extracts. It suits large uncompressed files and bundles of many internal files.
- **None mode** makes the server a catalog for the asset rather than a host. Users build it locally from its recipe.

## Choosing a mode

These are guidelines; use judgment for each asset and its users.

**Use file mode when:**

- The files are already compressed (`.gz`, `.bgz`, `.tbi`), so a tarball adds little.
- Each file is useful on its own, such as a `.bed.gz` someone wants without the rest.
- The asset has only a few files.
- The users are HTTP clients, such as WebAssembly frontends, cloud pipelines or REST APIs.

**Use archive mode when:**

- The asset is a large uncompressed file, such as a raw `.fa`, which compresses well.
- The asset has many files that only work as a set, such as STAR or bowtie2 index shards.
- The internal files are details of another tool rather than something users open.

**Use none mode when:**

- The server is a catalog of asset classes and recipes, not a data host.
- The asset is niche or experimental and users are expected to build it.
- You want to register the asset class before you host its data.

**Combine modes** (`[file, archive]`) when an asset is usually read file by file but is sometimes downloaded whole, for example browsed in a web UI and also pulled for local pipelines. The server then stores both forms, which doubles its disk use, so do this rarely.

### Seek key types as a hint

The seek key `type` in the asset class often points to the right mode:

- `type: file` means each seek key is a file users look up on its own. This fits **file mode**.
- `type: prefix` means a set of files that share a base name, such as the six `.bt2` files of a bowtie2 index. This fits **archive mode**.
- `type: directory` means the whole folder is the unit. This also fits **archive mode**.

## How modes affect other features

- **Remote hosting.** Any mode can be hosted on a remote server, CDN or S3 bucket. See [Configure remote storage](remotes.md).
- **DRS.** A file-mode asset is a [GA4GH DRS](https://ga4gh.github.io/data-repository-service-schemas/) bundle with one DRS object per file. An archive-mode asset is a single DRS object for the tarball.
- **Pull.** `refgenie pull` downloads and extracts the tarball for archive mode, downloads files one by one for file mode, and for none mode returns an error that explains how to build the asset locally.
- **Remote seek.** `refgenie seekr` can return a URL to one file of a file-mode asset without downloading it. Archive-mode assets must be pulled first.
