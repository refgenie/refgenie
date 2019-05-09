# <img src="img/refgenie_logo.svg" class="img-header"> a reference genome indexer

[![PEP compatible](http://pepkit.github.io/img/PEP-compatible-green.svg)](http://pepkit.github.io)

## What is refgenie?

Refgenie creates a standardized folder structure for reference genome files and indexes. You can download pre-built genomes or use the script to build your own for any genome you like. Minimal input is just a fasta file. It also comes with an optional [docker container (nsheff/refgenie)](https://hub.docker.com/r/nsheff/refgenie/) so you don't have to install any software.

## What makes refgenie better?

NGS pipelines require reference genomes for alignments and other computation. Usually, each pipeline has a unique way to organize reference genomes, and they may require the assembly to be organized in this structure. The problem with that is you must build a new reference genome folder organization for each pipeline, or at the very least specify your organization to each pipeline specifically. This makes it hard to share pipelines across environments and people, as they have different methods for organizing and passing index files to the pipelines.

If everyone built pipelines to follow a standard structure for reference genomes, then pipelines could just take a string describing that genome (e.g. "hg38") and would be able to know how to find the indexes it needs for that genome.

Refgenie standardizes the reference genome index folder structure, so you can build your pipelines around that standardized reference genome format. This makes it easy to switch to a new reference genome, or to share pipelines and reference data among collaborators. Pipelines then need only a single variable: the genome name, and you can even use an environment variable to make all pipelines work seamlessly. The goal is something like the [iGenomes project](http://support.illumina.com/sequencing/sequencing_software/igenome.html), except they just produce finalized files for you to download, but provide no software to produce your own reference for your own genomes. So, you can't use that to make a standardized format for your internal spike-in genomes or other species they don't provide.

You can download pre-computed tarballs of refgenie assemblies if you like, but more importantly, `refgenie` is a script, so you can produce the standard for whatever genome you want.

 

## Installing

Install from [GitHub releases](https://github.com/databio/refgenie/releases) or from PyPI with `pip`:


```console
pip install --user refgenie
```

Update with:

```console
pip install --user --upgrade refgenie
```

If the `looper` executable in not automatically in your `$PATH`, add the following line to your `.bashrc` or `.profile`:

```console
export PATH=~/.local/bin:$PATH
```

## Quick start

To test `looper`, follow the [Hello Looper example repository](https://github.com/databio/hello_looper) to run your first looper project:


```console
# download and unzip the hello_looper repository
wget https://github.com/databio/hello_looper/archive/master.zip
unzip master.zip

# Run looper:
cd hello_looper-master
looper run project/project_config.yaml
```

Detailed explanation of results is in the [Hello world tutorial](hello-world.md).
