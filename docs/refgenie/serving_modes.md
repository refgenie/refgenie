<h1>Serving modes</h1>

[TOC]

!!! tip "See also"
    To see serving modes in practice — including how to stage file-mode and archive-mode assets and serve them locally — see the [Building and Serving Tutorial](building_tutorial.md).

# Introduction

When you define an asset class, you specify what the asset *contains* through seek keys. But there is a separate question: how should a server *deliver* that asset to consumers? The `serving_modes` field on an asset class answers this question.

Different assets have fundamentally different access patterns. A FASTA file might be fetched individually by a browser-based tool, while a bowtie2 index -- a bundle of six tightly coupled shard files -- is only meaningful as a complete download. Serving modes let asset class authors express these access patterns so that servers can store and deliver data appropriately.

# The three serving modes

Refgenie defines three serving modes, each representing a different storage and delivery strategy:

| Mode | Server stores | What happens on request |
|------|--------------|------------------------|
| `file` | Extracted directory of individual files | Serves each file at its own URL |
| `archive` | `.tgz` tarball | Serves tarball for bulk download |
| `none` | Nothing (metadata only) | Returns metadata, digests, build commands |

An asset class specifies one or more of these as a list:

```yaml
serving_modes:
  - file
```

Or, less commonly:

```yaml
serving_modes:
  - file
  - archive
```

## File mode

In file mode, the server stores the asset as an extracted directory and exposes each file at its own URL. Consumers can fetch individual files with standard HTTP requests without downloading the entire asset.

This is the default mode and reflects the ecosystem's direction toward file-level API access as the primary pattern. It works well with HTTP-native consumers like web applications, cloud pipelines, and REST API clients.

## Archive mode

In archive mode, the server stores and serves the asset as a compressed `.tgz` tarball. Consumers download the entire archive and extract it locally.

This is the traditional refgenie distribution model. It suits assets that are large, uncompressed, or composed of many internal files that have no meaning in isolation.

## None mode

In none mode, the server does not host the asset data at all. It stores only metadata: the asset class definition, digests, and build commands. Consumers who need the data must build it locally using the associated recipe.

This mode turns the server into a catalog or registry rather than a data host.

# Guidelines for choosing a serving mode

The right serving mode depends on the characteristics of the asset. The following are guidelines, not hard rules -- use judgment based on the specific asset and its consumers.

## When to favor file mode

File mode is a good fit when:

- **Files are already compressed.** Formats like `.gz`, `.bgz`, and `.tbi` are already compressed individually. Wrapping them in a tarball adds packaging overhead without meaningful compression benefit.

- **Files are independently meaningful.** If consumers commonly want a specific file from the asset -- for example, fetching the `.bed.gz` without needing the rest -- file mode lets them do that with a single HTTP request.

- **The asset has a small number of files.** When an asset contains a handful of files, the overhead of individual HTTP requests is negligible compared to the flexibility gained.

- **Consumers are HTTP clients.** WebAssembly frontends, cloud pipelines, and REST APIs work naturally with URL-addressable files. File mode lets these consumers fetch exactly what they need without downloading and extracting archives.

## When to favor archive mode

Archive mode is a good fit when:

- **The asset is a large uncompressed file.** Raw FASTA files (`.fa`) are the canonical example -- they compress well, and tarball compression provides a real benefit for transfer.

- **The asset contains many files that are only meaningful as a complete set.** Aligner indices like STAR or bowtie2 produce dozens of internal shard files. Nobody ever wants one shard file from a STAR index; the files are opaque implementation details of the aligner.

- **Files are opaque to consumers.** When the internal file structure is an implementation detail rather than a user-facing interface, there is no benefit to individual file access.

## When to favor none mode

None mode is a good fit when:

- **The server is a catalog or registry** that documents available asset classes and their build recipes without hosting data.

- **The asset is niche or experimental** and users are expected to build it locally.

- **You want to register the asset class** and its metadata without committing to hosting data yet.

## Combined modes

An asset class can list multiple serving modes, such as `[file, archive]`. This tells the server to store the data in both forms, serving files individually for direct access and also providing a tarball for bulk download.

This comes at the cost of double disk usage on the server, so it should be uncommon. Most assets clearly fit one pattern. Consider combining modes when the asset is accessed file-by-file most of the time but occasionally needs bulk download -- for example, an asset that is both browsed in a web UI and pulled in bulk for local analysis pipelines.

# Seek key types as a signal

The seek key `type` field in the asset class definition is a useful indicator for which serving mode fits best:

- **`type: file`** seek keys indicate that the asset has individually addressable files. Each seek key points to a specific file that consumers look up independently. This maps naturally to **file mode**.

- **`type: prefix`** seek keys indicate that the asset is a set of sharded files sharing a common prefix. Bowtie2 indices are a typical example: six `.bt2` files that share a base name. The files are meaningless individually. This maps naturally to **archive mode**.

- **`type: directory`** seek keys indicate that the asset is an opaque directory bundle. The entire directory is the unit of interest, not any single file within it. This also maps naturally to **archive mode**.

This is not a mechanical rule. The seek key type reflects how the asset class author thinks about the asset's structure, and that same thinking informs the serving mode.

# How serving modes interact with other features

Serving modes are a property of the asset class, but they affect several parts of the refgenie ecosystem.

## Remote hosting

Serving modes are orthogonal to remote hosting. Any mode can be hosted on a remote server, CDN, or S3 bucket. The remote system handles URL redirection regardless of whether the underlying data is individual files or a tarball.

## DRS (GA4GH Data Repository Service)

File-mode assets map to DRS bundles with individually addressable child DRS objects -- each file becomes its own DRS object with its own identifier. Archive-mode assets map to a single DRS object representing the tarball.

## Pull

The `refgenie pull` command adapts its behavior to the serving mode:

- For **archive** mode, pull downloads the `.tgz` tarball and extracts it locally.
- For **file** mode, pull downloads files individually.
- For **none** mode, pull returns an error with guidance on how to build the asset locally.

## Seekr (remote seek)

Remote seek -- looking up an asset file by URL without downloading it first -- works with file-mode assets. The server can return a URL directly to the specific file. Archive-mode assets must be pulled and extracted before their files can be sought locally.

# Default

When an asset class does not specify `serving_modes`, it defaults to `[file]`. This reflects the ecosystem's direction toward file-level API access as the primary consumption pattern for reference genome assets.
