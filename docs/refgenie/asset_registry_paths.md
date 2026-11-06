# Asset registry paths

An *asset registry path* names one asset, or one file inside it, on the command line:

```console
genome/asset_group.seek_key:asset
```

| Part | Meaning | If you leave it out |
|---|---|---|
| `genome` | A genome alias, such as `hg38` | Required (or pass `--genome-digest` instead) |
| `asset_group` | The group of related assets, such as `fasta` or `bowtie2_index` | Required |
| `seek_key` | One file or value inside the asset, such as `fai` | The asset's default seek key |
| `asset` | The asset's name within the group, such as `default` or `2.5.4` | The group's default asset |

The default seek key is the asset's only path seek key, or, if it has several, the one with the same name as the asset group. The group's default asset is the one you built most recently, or, if you only pulled, the first one you pulled.

These commands all return the same path:

```console
$ refgenie seek rCRSd/fasta
/path/to/genomes/alias/rCRSd/fasta/default/rCRSd.fa

$ refgenie seek rCRSd/fasta.fasta
/path/to/genomes/alias/rCRSd/fasta/default/rCRSd.fa

$ refgenie seek rCRSd/fasta.fasta:default
/path/to/genomes/alias/rCRSd/fasta/default/rCRSd.fa
```

In files that `refgenie populate` reads, write the path with a `refgenie://` prefix, as in `refgenie://hg38/fasta`. See [Populate](populate.md).
