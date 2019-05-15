# <img src="img/refgenie_logo.svg" class="img-header"> genome index manager

[![PEP compatible](http://pepkit.github.io/img/PEP-compatible-green.svg)](http://pepkit.github.io)

## What is refgenie?

Refgenie creates a standardized folder structure for reference genome files and indexes. You can download pre-built genomes or build your own for any fasta file.

## What makes refgenie better?

Refgenie provides a **standard folder structure** for reference genome indexes, so that alignment tools can easily swap from one genome to another. Most importantly, Refgenie is **scripted** so that users can create their own assembly index packages from whatever genome source they like.

## Installing

If you just want to use pre-built refgenie assemblies, just head over to the [download page](download.md); you don't even need to install refgenie. If you want to index your own genomes, then you'll need to install refgenie plus your genome indexers of choice. Install refgenie from [GitHub releases](https://github.com/databio/refgenie/releases) or from PyPI with `pip`:


```console
pip install --user refgenie
```

Update with:

```console
pip install --user --upgrade refgenie
```

After that, you'll need to [install the genome indexers](install.md) -- but first, you can confirm that `refgenie` is functioning:

## Quick start

See if your install worked by invoking `refgenie` from the command line:

```
refgenie -h
```

If the `refgenie` executable in not automatically in your `$PATH`, add the following line to your `.bashrc` or `.profile` (or `.bash_profile` on MACOS):

```console
export PATH=~/.local/bin:$PATH
```

Next, [install the genome indexers](install.md).