# <img src="img/refgenie_logo.svg" class="img-header"> a reference genome indexer

[![PEP compatible](http://pepkit.github.io/img/PEP-compatible-green.svg)](http://pepkit.github.io)

## What is refgenie?

Refgenie creates a standardized folder structure for reference genome files and indexes. You can download pre-built genomes or use the script to build your own for any genome you like. Minimal input is just a fasta file. It also comes with an optional [docker container (nsheff/refgenie)](https://hub.docker.com/r/nsheff/refgenie/) so you don't have to install any software.

## What makes refgenie better?

Refgenie provides a **standard folder structure** for reference genome indexes, so that alignment tools can easily swap from one genome to another. Most importantly, Refgenie is **scripted** so that users can create their own assembly index packages from whatever genome source they like.

## Installing

Install from [GitHub releases](https://github.com/databio/refgenie/releases) or from PyPI with `pip`:


```console
pip install --user refgenie
```

Update with:

```console
pip install --user --upgrade refgenie
```

If the `refgenie` executable in not automatically in your `$PATH`, add the following line to your `.bashrc` or `.profile`:

```console
export PATH=~/.local/bin:$PATH
```

## Quick start

To test `refgenie`, follow the [Hello Looper example repository](https://github.com/databio/hello_looper) to run your first looper project:


```console
# download and unzip the hello_looper repository
wget https://github.com/databio/hello_looper/archive/master.zip
unzip master.zip

# Run looper:
cd hello_looper-master
looper run project/project_config.yaml
```

Detailed explanation of results is in the [Hello world tutorial](hello-world.md).
