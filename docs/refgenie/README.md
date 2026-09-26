
<div class="page-hero">
<img src="img/refgenie_logo.svg" alt="Refgenie" class="page-hero-logo" />
<p class="page-hero-tagline">reference genome manager</p>
<div class="badges">
<a href="https://pepkit.github.io/img/PEP-compatible-green.svg"><img src="https://pepkit.github.io/img/PEP-compatible-green.svg"></a>
</div>
</div>


## What is refgenie?

Refgenie manages storage, access, and transfer of reference genome resources. It provides command-line and Python interfaces to *download* pre-built reference genome "assets", like indexes used by bioinformatics tools. It can also *build* assets for custom genome assemblies. Refgenie provides programmatic access to a standard genome folder structure, so software can swap from one genome to another.

## What makes refgenie better?

1. **It provides a command-line interface to download individual resources**. Think of it as `GitHub` for reference genomes. You just type `refgenie pull hg38/bwa_index`.

2. **It's scripted**. In case you need resources *not* on the server, such as for a custom genome, you can `build` your own: `refgenie build custom_genome/bowtie2_index`.

3. **It simplifies finding local asset locations**. When you need a path to an asset, you can `seek` it, making your pipelines portable across computing environments: `refgenie seek hg38/salmon_index`.

4. **It provides remote operation mode**, useful for cloud applications. Get a path to an asset file hosted on AWS S3: `refgenie seekr hg38/fasta --remote-class s3`.

5. **It includes a Python API**. For tool developers, you use `from refgenie import Refgenie` to get a Python object with paths to any genome asset, *e.g.*, `rgc = Refgenie(); rgc.asset.seek("hg38", "kallisto_index")`.

6. **It strictly determines genomes compatibility**. Users refer to genomes with arbitrary aliases, like "hg38", but refgenie uses sequence-derived identifiers to verify genome identity with asset servers.

7. **It is scalable for large-scale operations**. Refgenie is backed by a database, enabling efficient, multi-user, management of genome assets. By default, it uses a local SQLite database, but it can be configured to use PostgreSQL, including setups on remote servers, to support high-performance and distributed workflows.

## Quick example

### Install

Refgenie is a Python package, install from [PyPi](https://pypi.org/project/refgenie/):

```console
pip install refgenie
```

And that's it! If you wish to use refgenie in *remote mode*, see [further reading on remote mode in refgenie](remotes.md).

If you're connected to the Internet, call a test command, e.g.:

```console
refgenie seekr hg38/fasta
```

### Initialize to use refgenie locally

By default, Refgenie keeps track of what's available using local configuration initialized by `refgenie init`:

```console
refgenie init
```

See [further reading on configuring refgenie](configuration.md).

### Download indexes and assets for a remote reference genome

Use `refgenie pull` to download pre-built assets from a remote server. View available remote assets with `listr`:

```console
refgenie listr
```

Response:
```console
                 Refgenie assets. Source: http://refgenomes.databio.org
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━┓
┃ Genome digest                                    ┃ Asset group             ┃ Asset   ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━┩
│ 047c6e1eda552b50c5add59ff0995a40bc4ce1732e3cc4ae │ bowtie2_index           │ default │
│ 047c6e1eda552b50c5add59ff0995a40bc4ce1732e3cc4ae │ bwa_index               │ default │
│ 047c6e1eda552b50c5add59ff0995a40bc4ce1732e3cc4ae │ fasta                   │ default │
....
```

Next, pull one:

```console
refgenie pull rCRSd/bowtie2_index
```

Response:
```console
WARNING  No local digest for genome alias: rCRSd. Setting genome identity with server: http://refgenomes.databio.org
INFO     Connected to server: title='refgenieserver' version='0.7.0'
INFO     Setting 'rCRSd' identity with server: http://refgenomes.databio.org
INFO     Determined digest for rCRSd: 94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabf8c1c2db5a6c0c7b2a8
INFO     Set genome alias: rCRSd
rCRSd/bowtie2_index:default ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100% • 52.3/52.3 kB
...
```

See [further reading on downloading assets](pull.md).

### Build your own indexes and assets for a custom reference genome

Refgenie assets are scripted, so if what you need is not available remotely, you can use `build` it locally:


```console
refgenie build mygenome/bwa_index
```

See [further reading on building assets](build.md).

### Retrieve paths to *local* refgenie-managed assets

Once you've populated your refgenie with a few assets, use `seek` to retrieve their local file paths:

```console
refgenie seek mm10/bowtie2_index
```

This will return the path to the particular asset of interest, regardless of your computing environment. This gives you an ultra-portable asset manager! See [further reading on retrieving asset paths](seek.md).

### Retrieve paths to *remote* refgenie-managed assets

Use `seekr` (short for "seek remote") to retrieve remote `seek_key` targets:

```console
refgenie seekr mm10/fasta.fai
```

This will return the path to the particular remote file of interest, here: FASTA index file, which is a part of `mm10/fasta` asset.

See [further reading on using refgenie in remote mode](remotes.md).

### Use refgenie from Python

```python
from refgenie import Refgenie

rgc = Refgenie()
rgc.asset.seek("hg38", "bowtie2_index")
```

The `Refgenie` object connects to the configured database and provides the same functionality as the CLI. See [the Refgenie Python object](refgenie_object.md) for details.

---
If you want to read more about the motivation behind refgenie and the software engineering that makes refgenie work, proceed next to the [overview](overview.md).
