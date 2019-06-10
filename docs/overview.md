# Refgenie overview

## Motivation

Enormous effort goes into assembling and curating reference genomes. Reference assemblies are the starting point for many downstream tools, such as for sequence alignment and annotation. Many tools produce independent assets that accompany a genome assembly; for instance, many aligners like `bwa` or `bowtie2` must *hash* the genome, creating *indexes* to improve alignment performance. These indexes are typically shared among many pipelines and tools, so it's common practice for a research group to organize a central folder for reference genome assets, which includes indexes and other files like annotations. In addition to saving disk space, this simplifies sharing software among group members. However, each group typically does this independently.

<img src="../img/refgenie_interfaces.svg" style="float:right; width:350px">

One effort to distribute standard, organized reference sequences and annotation files is Illumina's *IGenomes* project, which distributes pre-built archives of common assets for common genomes. But IGenomes doesn't provide automated software to produce a standard reference for an arbitrary genome. This is problematic because we often need to align to a custom genome, such as a spike-in control or a personal genome assembly. Furthermore, packaging many resources together in a single archive disallows modular access to individual genome assets, costing computational resources.

**Refgenie simultaneously provides structure to manually build assets while improving modular access to pre-built assets in the same system. Refgenie does this by providing two ways to obtain genome assets** (see figure right): first, web-based access to individual pre-built assets via web interface or application programming interface (API); and second, an interface for scripted asset builds that reproduce structured resource output for arbitrary genome inputs. This enables users to either retrieve or produce identically structured outputs on demand for any genome of interest, including new assemblies, private assemblies, or custom genomes for which a public set of assets cannot exist. 

## Refgenie ecosystem

Refgenie consists of 3 independent tools that work together to enable the above scenario:

### 1. The `refgenie` command-line interface (CLI)

The component that most users will interact with is the command-line interface. With a simple `pip install`, this provides the `refgenie` command, which can be used to `pull` or `build` an asset of interest.

### 2. The `refgenieserver` package.

While the `refgenie` CLI has utility along, and a user would achieve great benefit just by using `refgenie build` with no remote component, it becomes even more powerful when it can communite with a remote server via `refgenie pull`. In order to do this, we require a remote counterpart to fill the `pull` requests, and that is handled by the `refgenieserver` software. Refgenieserver provides a web interface and an API that can be used by the CLI (or by any other tool as well), which allows users to list available assets and download them.

### 3. The `refgenconf` python package

Both the server and the CLI rely on a separate package that has additional utility: the `refgenconf` python package. This package provides a python API to interact with the genome configuration files produced by refgenie. This provides a simple interface whereby any third-party python packages can have easy access to the refgenie organizational structure.
