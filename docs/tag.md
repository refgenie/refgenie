# Tagging assets

## Why to tag assets?
It is natural in a research environment to use various flavors of the reference genome related resources that may result from different versions of the software used to create them. And this is what inspired the introduction of assets tagging concept in `refgenie`.

Tag can be **any text or number**, so it is well suited to contain software version information or even a concise description, like `0.4.1` or `new_build_strategy`

## How to tag assets?

Asset tagging was designed to be very flexible -- assets can be tagged at any point:

- **tagging when assets are built:**

```console
export REFGENIE="genome_config.yaml"
refgenie init -c $REFGENIE 
refgenie pull hg38/fasta
refgenie build hg38/bowtie2_index:2.3.5.1
```
or
```console
refgenie build hg38/bowtie2_index:2.3.3.1
```

- **tagging already built/pulled assets (re-tagging):**

```console
refgenie tag hg38/bowtie2_index:2.3.5.1 --tag most_recent
```

- **no tagging at all:**

Importantly, you don't have to care about tags at all if you don't need to because there is a **default** representative for every asset in your assets inventory.

Building without specifying a tag will tag the asset as `default`

```console
refgenie build hg38/bwa_index
```

### Default tags

If you pull or build the first asset of a given kind it will become the default one, which `refgenie` will use for any actions when no tag is explicitly specified. For example the

```console
refgenie seek hg38/bowtie2_index 
```

call would return the path to the asset tagged with `most_recent` since it was the first `bowtie2_index` asset built/pulled for `hg38` genome.

To retrieve the path to any other asset, you need to specify the tag:

```console
refgenie seek hg38/bowtie2_index:2.3.3.1 
``` 