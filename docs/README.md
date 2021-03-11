
# <img src="img/refgenie_logo.svg" class="img-header"> reference genome manager

[![PEP compatible](https://pepkit.github.io/img/PEP-compatible-green.svg)](https://pepkit.github.io)
[![PyPi](https://img.shields.io/pypi/v/refgenie.svg)](https://pypi.org/project/refgenie/)

<iframe src="https://ghbtns.com/github-btn.html?user=refgenie&repo=refgenie&type=star&count=true" frameborder="0" scrolling="0" width="150" height="20" title="GitHub"></iframe>

## What is refgenie?

Refgenie manages storage, access, and transfer of reference genome resources. It provides command-line and Python interfaces to *download* pre-built reference genome "assets", like indexes used by bioinformatics tools. It can also *build* assets for custom genome assemblies. Refgenie provides programmatic access to a standard genome folder structure, so software can swap from one genome to another.

## What makes refgenie better?

1. **It provides a command-line interface to download individual resources**. Think of it as `GitHub` for reference genomes. You just type `refgenie pull hg38/bwa_index`.

2. **It's scripted**. In case you need resources *not* on the server, such as for a custom genome, you can `build` your own: `refgenie build custom_genome/bowtie2_index`.

3. **It simplifies finding local asset locations**. When you need a path to an asset, you can `seek` it, making your pipelines portable across computing environments: `refgenie seek hg38/salmon_index`.

4. **It includes a python API**. For tool developers, you use `rgc = refgenconf.RefGenConf("genomes.yaml")` to get a Python object with paths to any genome asset, *e.g.*, `rgc.seek("hg38", "kallisto_index")`.

5. **It strictly determines genomes compatibility**. Users refer to genomes with arbitrary aliases, like "hg38", but refgenie uses sequence-derived identifiers to verify genome identity with asset servers.


## Quick example

### Install and initialize

Refgenie keeps track of what's available using a configuration file initialized by `refgenie init`:

```console
pip install --user refgenie
export REFGENIE='genome_config.yaml'
refgenie init -c $REFGENIE
```

### Download indexes and assets for a remote reference genome

Use `refgenie pull` to download pre-built assets from a remote server. View available remote assets with `listr`:

```console
refgenie listr
```

Response:
```console
                        Remote refgenie assets                        
                 Server URL: http://refgenomes.databio.org                 
┏━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ genome              ┃ assets                                       ┃
┡━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ mouse_chrM2x        │ fasta, bwa_index, bowtie2_index              │
│ hg38                │ fasta, bowtie2_index                         │
│ rCRSd               │ fasta, bowtie2_index                         │
│ human_repeats       │ fasta, hisat2_index, bwa_index               │
└─────────────────────┴──────────────────────────────────────────────┘
```

Next, pull one:

```console
refgenie pull rCRSd/bowtie2_index
```

Response:
```console
Downloading URL: http://rg.databio.org/v3/assets/archive/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/bowtie2_index
94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/bowtie2_index:default ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 100.0% • 128.0/117.0 KB • 1.8 MB/s • 0:00:00
Download complete: /Users/mstolarczyk/Desktop/testing/refgenie/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/bowtie2_index/bowtie2_index__default.tgz
Extracting asset tarball: /Users/mstolarczyk/Desktop/testing/refgenie/data/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/bowtie2_index/bowtie2_index__default.tgz
Default tag for '94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/bowtie2_index' set to: default
Created alias directories: 
 - /Users/mstolarczyk/Desktop/testing/refgenie/alias/rCRSd/bowtie2_index/default
```

See [further reading on downloading assets](pull.md).

### Build your own indexes and assets for a custom reference genome

Refgenie assets are scripted, so if what you need is not available remotely, you can use `build` it locally:


```console
refgenie build mygenome/bwa_index
```

See [further reading on building assets](build.md).

### Retrieve paths to refgenie-managed assets

Once you've populated your refgenie with a few assets, use `seek` to retrieve their local file paths:

```console
refgenie seek mm10/bowtie2_index
```

This will return the path to the particular asset of interest, regardless of your computing environment. This gives you an ultra-portable asset manager! See [further reading on retrieving asset paths](seek.md).

If you want to read more about the motivation behind refgenie and the software engineering that makes refgenie work, proceed next to the [overview](overview.md).
