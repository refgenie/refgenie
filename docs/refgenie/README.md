
<div class="page-hero">
<img src="img/refgenie_logo.svg" alt="Refgenie" class="page-hero-logo" />
<p class="page-hero-tagline">reference genome manager</p>
<div class="badges">
<a href="https://pepkit.github.io/img/PEP-compatible-green.svg"><img src="https://pepkit.github.io/img/PEP-compatible-green.svg"></a>
</div>
</div>


## What is refgenie?

Refgenie manages storage, access, and transfer of reference genome resources. It provides command-line and Python interfaces to *download* pre-built reference genome "assets", like indexes used by bioinformatics tools. It can also *build* assets for custom genome assemblies. Refgenie provides programmatic access to a standard genome folder structure, so software can swap from one genome to another.

!!! warning "Alpha"
    Refgenie 1.0 is an alpha. It works, but commands and data can still change, so pin the version in production pipelines (for example, `refgenie==1.0.0a2`).

## What makes refgenie better?

1. **It downloads individual resources from the command line**. Think of it as `GitHub` for reference genomes: `refgenie pull hg38/bwa_index`.

2. **It's scripted**. If you need resources that are not on a server, such as for a custom genome, you can `build` your own: `refgenie build custom_genome/bowtie2_index`.

3. **It finds local asset locations**. When you need a path to an asset, you `seek` it, so your pipelines are portable across computing environments: `refgenie seek hg38/salmon_index`.

4. **It works with remote files**, which is useful for cloud applications. Get the URL of an asset file on a server without downloading it: `refgenie seekr hg38/fasta`.

5. **It includes a Python API**. Tool developers use `from refgenie import Refgenie` to get paths to any genome asset from Python.

6. **It checks genome compatibility**. You refer to genomes by aliases like "hg38", but refgenie identifies each genome by a digest computed from its sequences ([GA4GH sequence collections](https://ga4gh.github.io/seqcol-spec/)), so assets always match the genome they were built from.

7. **It scales**. Refgenie keeps its records in a database. It uses a local SQLite file by default and can use PostgreSQL, including on a remote server, for multi-user and distributed setups.

## Quick start

Four commands install refgenie, set it up, download a small asset, and print its path:

```console
pip install "refgenie>=1.0.0a2"
refgenie init
refgenie pull hg38/fasta_index
refgenie seek hg38/fasta_index.chrom_sizes
```

1. **Install.** Version 1.0 is a pre-release, so ask [PyPI](https://pypi.org/project/refgenie/) for it by version. See [Install](install.md).
2. **Initialize.** `refgenie init` creates the local database and genome folder, and subscribes to the public server at `https://api.refgenie.org`. See [Configure](configuration.md) for other locations and database options.
3. **Pull.** `hg38/fasta_index` is about 12 KB: the FASTA index and chromosome sizes for hg38, without the 3 GB FASTA. Not sure which hg38 you need? See [Which hg38 should I use?](choosing_hg38.md).
4. **Seek.** `refgenie seek` prints the local path to a file in the asset, so your scripts work on any computer. See [Retrieve paths to assets](seek.md).

To see everything you can pull, run `refgenie listr`, or browse [refgenie.org](https://refgenie.org). See [Download pre-built assets](pull.md).

## Next steps

### Build assets for a custom genome

If what you need is not on a server, register your genome from its FASTA file and build assets locally. Building needs recipes from a [data channel](data_channels.md):

```console
refgenie data-channel add refgenie https https://refgenie.github.io/refgenie-registry/index.yaml
refgenie data-channel sync refgenie --exists-ok
refgenie genome init --fasta mygenome.fa.gz --name mygenome
refgenie build mygenome/bwa_index
```

See [Build assets](build.md).

### Get paths to remote assets

Use `seekr` ("seek remote") to get the URL of a file on a server, here the FASTA index inside the `hg38/fasta` asset:

```console
refgenie seekr hg38/fasta.fai
```

See [Configure remote storage](remotes.md).

### Use refgenie from Python

```python
from refgenie import Refgenie

rgc = Refgenie()
digest = rgc.alias.resolve("hg38")  # seek takes a genome digest, not an alias
rgc.asset.seek(digest, "fasta_index", seek_key_name="chrom_sizes")
```

The `Refgenie` object connects to the configured database and does everything the CLI does. See [the Refgenie Python object](refgenie_object.md).

---
To read more about why refgenie exists and how it is built, continue to the [overview](overview.md).
