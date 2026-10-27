# Tagging assets

## Why to tag assets?

It is natural in a research environment to use various flavors of the reference genome related resources that may result from different versions of the software used to create them. And this is what inspired the introduction of assets tagging concept in `refgenie`.

### Tag character whitelist

Tag can be **any text or number** composed of characters safe for Uniform Resource Identifiers (URIs) as per [RFC3986](https://www.ietf.org/rfc/rfc3986.txt), so it is well suited to contain software version information or even a concise description, like `0.4.1` or `new_build_strategy`.

RFC3986 section 2.3. _Unreserved Characters_: characters that are allowed in a URI include uppercase and lowercase letters, decimal digits, hyphen, period, underscore, and tilde.

```
"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-._~"
```

## How to tag assets?

Asset tagging is very flexible. You can tag assets when you build them, add or change tags to already built assets, or just not use tags at all if you don't need them.

### Tagging when assets are built

Here we'll demonstrate how you can specify a tag when building an asset:

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

### Tagging already built/pulled assets (re-tagging)

If you already built an asset, you can add a tag to it. Here, we'll add a tag for `most_recent` to our bowtie2 index asset:

```console
refgenie tag hg38/bowtie2_index:2.3.5.1 --tag most_recent
```

Now you could retrieve this asset using either of those tags. In other words, _assets can have more than 1 tag_.

### No tagging at all

Importantly, you don't have to care about tags at all if you don't need to because there is a **default** tag for every asset in your assets inventory. Building without specifying a tag will tag the asset as `default`. If you don't specify a tag when trying to retrieve an asset path, it will assume you're looking for the default tagged asset.

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

### Changing the default tag

If you want to make a tag the default one, use the `-d`/`--default` option in `refgenie tag` command:

```console
refgenie tag hg38/bowtie2_index:2.3.3.1 -d
```
