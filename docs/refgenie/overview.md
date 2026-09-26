# Refgenie overview

## Motivation

Reference genome assemblies are the starting point for many downstream tools, such as sequence alignment and annotation. Many tools produce independent assets that accompany a genome assembly; for instance, aligners like `bwa` or `bowtie2` must *hash* the genome, creating *indexes* to improve alignment performance. These indexes are typically shared across tools, so it's common for a research group to organize a central folder for reference genome assets, which includes indexes and other files like annotations. In addition to saving disk space, this simplifies sharing software among group members. Structuring assets uniformly allows software to adapt from one reference assembly to the next. However, each group typically does this independently. If we could standardize this across groups, it would make it easier to share scripts and software that use genome-related resources.

One effort to distribute standard, organized reference sequences and annotation files is Illumina's *IGenomes* project, which distributes pre-built archives of common assets for common genomes. This allows multiple groups to share one standard, but it has a few limitations: IGenomes doesn't provide software to produce a standard reference for an arbitrary genome. This is problematic because we often need to align to a custom genome, such as a spike-in control or a personal genome assembly. Furthermore, packaging many resources together in a single archive precludes itemized access to individual genome assets, costing computational resources.

<img src="../img/refgenie_interfaces.svg" style="float:right; width:350px">

### Functionality

 **Refgenie simultaneously provides structure to manually build assets while improving modular access to pre-built assets in the same system.** Refgenie does this by providing two ways to obtain genome assets (see figure at right).

  1. Web-based access to individual pre-built assets via web interface or application programming interface (API)
  2. An interface for scripted asset "builds," each of which produces structured output for arbitrary genome inputs.

This two-pronged approach enables users to either retrieve or produce *identically structured* outputs on demand for *any genome* of interest, including new assemblies, private assemblies, or custom genomes for which a public set of assets cannot exist.

## Refgenie ecosystem

Refgenie is a single unified package that provides everything you need to manage, build, serve, and distribute reference genome assets.

### The `refgenie` command-line interface (CLI)

A simple `pip install refgenie` provides the `refgenie` command, which can be used to `pull` or `build` an asset of interest. Additionally, `refgenie` can be used programmatically from Python through its Python API, allowing developers to integrate refgenie functionality directly into their scripts and applications. For detailed usage examples and API documentation, see the [refgenie usage documentation](code/refgenie.md).

### Built-in server

Refgenie includes a built-in server that you start with `refgenie serve`. The server provides a web interface and REST API that can be used by the CLI (or by any other tool), allowing users to list available assets and download them. We host a public instance at [refgenomes.databio.org](http://refgenomes.databio.org), but you can also run your own instance. See the [server documentation](server/README.md) for details.

### Two operating modes

Refgenie supports two database backends:

| Backend | Use case | Setup |
|---|---|---|
| **SQLite** (default) | Single-user workstations, local development | Zero configuration -- works out of the box |
| **PostgreSQL** | Multi-user servers, production deployments, cloud | Requires a PostgreSQL server; configured via a YAML file |

See [Database backends](database.md) for configuration details.

### RefgetStore

When you initialize a genome with `refgenie genome init`, refgenie loads all sequences into a local **RefgetStore** -- a content-addressable store where every sequence is identified by its GA4GH digest. This enables deduplication across genomes, direct subsequence retrieval with `refgenie getseq`, and verification of genome identity. See [Genome initialization and the RefgetStore](genome_tutorial.md) for the full workflow.

### Data channels

Refgenie's extensibility is powered by **data channels** -- external sources that provide definitions for asset classes and recipes. Data channels enable the community to publish, distribute, and synchronize new types of reference genome assets and build instructions without requiring changes to the core refgenie codebase. See [Data channels](data_channels.md) for details.

### Staging and remote storage

Built assets can be **staged** for serving and then **pushed** to remote storage like S3. Staging prepares assets according to their serving modes (file or archive). See [Stage assets for serving](staging.md) and [Configure remote storage](remotes.md).
