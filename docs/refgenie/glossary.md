# Glossary

- **Genome** (also *assembly*, *reference genome*): a set of sequences for an organism, such as hg38.
- **Genome digest**: a genome's identifier, the GA4GH [sequence collection](https://ga4gh.github.io/seqcol-spec/) digest of its sequences. For example, `jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP`.
- **Alias**: a human-readable name for a genome digest, such as `hg38`. A genome can have several aliases.
- **Asset group**: a kind of asset for a genome, such as `fasta` or `bowtie2_index`. A group holds one or more assets, and one is the default.
- **Asset**: one named build within an asset group, such as `default` or `2.5.5`, made of one or more files.
- **Asset class**: the definition of a kind of asset: its files and seek keys.
- **Recipe**: the instructions for building an asset of a given class.
- **Data channel**: a published collection of asset classes and recipes that you can sync.
- **Seek key**: the name of one file, folder, or value inside an asset, such as `fai` in the `fasta` asset. `refgenie seek` uses it to return a path.
- **Asset registry path**: a string that names an asset, of the form `genome/asset_group.seek_key:asset`, such as `hg38/fasta.fai:default`. The seek key and asset are optional. See [Refer to assets](asset_registry_paths.md).
- **RefgetStore**: the local store of sequences, each saved under its GA4GH digest.
- **Server**: a refgenie instance that serves assets over HTTP, such as `https://api.refgenie.org`.
- **Remote**: a storage location, such as an S3 bucket, that staged assets are pushed to.
