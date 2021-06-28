# Glossary

A common set of concepts is central to these and other genomics resources docs,
but often terminology differs between sources, or even within the same one (and may here).

Here are some clarifications regarding wording and meaning.

- **Assembly**: a particular version of a consensus genomic sequence for an organism
- **Asset**: a folder consisting of one or more files related to a specific reference genome
- **Asset registry path**: a string used to refer to assets, of the form: `genome/asset:tag`. If using an asset with more than 1 seek key, additional keys are appended like: `genome/asset:key:tag`.
- **Genome**: here, this is used interchangeably with *assembly*
- **Reference assembly**: another synonym for *assembly*
- **Reference genome**: yet another synonym for assembly
- **Seek key** a string identifier for a particular file or folder contained within an asset. A seek key is used to retrieve a path with refgenie seek. Seek keys allow a single asset to provide more than one endpoint for retrieval.
- **Tag**: a unique string identifier that allows for multiple assets of the same name to co-exist. One common use case for tags is to maintain multiple versions of an asset.
