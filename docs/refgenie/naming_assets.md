# Asset groups and naming

An *asset group* holds several assets of the same kind for one genome, such as `bowtie2_index` assets built with different bowtie2 versions or settings. Each asset in the group has its own *name*.

## Asset names

Names can use letters, digits, hyphen, period, underscore, and tilde (the URI-safe characters of [RFC 3986](https://www.rfc-editor.org/rfc/rfc3986#section-2.3)), so they can hold a version or a short label, like `2.5.4` or `new_build_strategy`.

Most recipes name an asset after the version of the tool that built it. To choose the name yourself, add it after a colon when you build:

```console
refgenie build hg38/bowtie2_index:my_custom_name
```

To rename an asset you already have:

```console
refgenie rename hg38/bowtie2_index:my_custom_name -n most_recent
```

## The default asset

Each group has a *default* asset, which refgenie uses when you leave out the name. The asset you build becomes the default. A pulled asset becomes the default only if it is the first one in its group.

```console
refgenie seek hg38/bowtie2_index                  # the default asset
refgenie seek hg38/bowtie2_index:my_custom_name   # a named asset
```

## More than one group of a kind

A genome can also hold several groups of the same [asset class](asset_class_specification.md), such as a `fasta` group for the genome and a `fasta_txome` group for its transcriptome. Each recipe input names the group it uses by default (`bowtie2_index` uses `fasta`; `salmon_index` uses `fasta_txome`). Use `--assets` to pick another:

```console
refgenie build hg38/bowtie2_index:txome --assets fasta=hg38/fasta_txome
```
