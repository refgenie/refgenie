# Managing asset groups and naming

## What are asset groups?

Asset groups in refgenie allow you to organize multiple related assets of the same type within a single group. This is useful when you have different versions or variants of the same asset type that result from different software versions, build parameters, or methodologies.

For example, you might have multiple `bowtie2_index` assets built with different versions of bowtie2, or multiple `fasta` assets representing different assemblies or filtered versions of the same genome.

## Why use asset groups?

It is natural in a research environment to use various flavors of reference genome-related resources that may result from different versions of the software used to create them, or different processing approaches. Asset groups provide a clean way to organize these variants while maintaining clear identification and easy access.

Benefits of asset groups include:

- **Version management**: Keep multiple versions of the same asset type
- **Method comparison**: Store assets built with different tools or parameters
- **Flexibility**: Switch between different variants as needed
- **Organization**: Maintain a clean, organized asset structure

## Asset naming

### Asset name character requirements

Asset names can contain **any text or number** composed of characters safe for Uniform Resource Identifiers (URIs) as per [RFC3986](https://www.ietf.org/rfc/rfc3986.txt), making them well suited to contain software version information or concise descriptions, like `v2.3.5.1` or `new_build_strategy`.

RFC3986 section 2.3. _Unreserved Characters_: characters that are allowed in a URI include uppercase and lowercase letters, decimal digits, hyphen, period, underscore, and tilde.

```console
"abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-._~"
```

### Naming assets when building

When building assets, you can specify a custom name to organize them within asset groups:

```console
refgenie init
refgenie pull hg38/fasta
refgenie build hg38/bowtie2_index:my_custom_name
```

### Renaming existing assets

You can rename assets that have already been built or pulled using the `rename` command:

```console
refgenie rename hg38/bowtie2_index -n most_recent
```

This allows you to update asset names to reflect their purpose or version after they've been created.

### Default asset behavior

Every asset group has a **default** asset that refgenie will use when no specific asset name is provided. When you build or pull the first asset of a given type, it becomes the default asset for that group.

```console
refgenie build hg38/bwa_index
```

This creates a default `bwa_index` asset.

### Accessing specific assets

To access the default asset in a group, simply omit the asset name:

```console
refgenie seek hg38/bowtie2_index
```

To access a specific named asset within the group, specify the name:

```console
refgenie seek hg38/bowtie2_index:my_custom_name
```
