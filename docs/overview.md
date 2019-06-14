# Refgenie overview

## Motivation

Enormous effort goes into assembling and curating reference genomes. Reference assemblies are the starting point for many downstream tools, such as for sequence alignment and annotation. Many tools produce independent assets that accompany a genome assembly; for instance, many aligners like `bwa` or `bowtie2` must *hash* the genome, creating *indexes* to improve alignment performance. These indexes are typically shared among many pipelines and tools, so it's common practice for a research group to organize a central folder for reference genome assets, which includes indexes and other files like annotations. In addition to saving disk space, this simplifies sharing software among group members. However, each group typically does this independently.

<img src="../img/refgenie_interfaces.svg" style="float:right; width:350px">

One effort to distribute standard, organized reference sequences and annotation files is Illumina's *IGenomes* project, which distributes pre-built archives of common assets for common genomes. But IGenomes doesn't provide automated software to produce a standard reference for an arbitrary genome. This is problematic because we often need to align to a custom genome, such as a spike-in control or a personal genome assembly. Furthermore, packaging many resources together in a single archive precludes itemized access to individual genome assets, costing computational resources.

### Functionality

Uniformly structuring assets regardless of reference assembly allows both standard and custom software to be used and/or written without regard for idiosyncrasies of a particular context. **Refgenie simultaneously provides *structure to manually build assets* while improving modular access to pre-built assets in the same system. Refgenie does this by providing *two ways to obtain genome assets*** (see figure at right).
  1. Web-based access to individual pre-built assets via web interface or application programming interface (API)
  2. An interface for scripted asset "builds," each of which produces structured output for arbitrary genome inputs. This enables users to either retrieve or produce *identically structured* outputs on demand for *any genome* of interest, including new assemblies, private assemblies, or custom genomes for which a public set of assets cannot exist. 

## Refgenie ecosystem

Refgenie consists of 3 independent tools that work together to provide the [intended functionality](#functionality):

### 1. The `refgenie` command-line interface (CLI)

The component that most users will interact with is the command-line interface. A simple `pip install refgenie` provides the `refgenie` command, which can be used to `pull` or `build` an asset of interest.

### 2. The `refgenieserver` package

While the `refgenie` CLI is useful by itself, and there's great benefit with simply `refgenie build` without any remote component, it becomes even more powerful when it can communicate with a remote server via `refgenie pull`. In order to do this, we require a remote counterpart to fulfill the `pull` requests. That is handled by `refgenieserver`, which provides a web interface and an API that can be used by the CLI (or by any other tool), and allows users to list available assets and download them.

### 3. The `refgenconf` Python package

Both the server and the CLI rely on `refgenconf`. This package provides a Python API to interact with the genome configuration files produced by Refgenie. This provides a simple interface whereby any third-party Python packages can easily leverage the asset organizatoion that Refgenie implements.
