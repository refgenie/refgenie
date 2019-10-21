# Refgenie genome configuration file

Refgenie will read and write a genome configuration file in yaml format. In general, you shouldn't need to mess with the config file. You create one with `refgenie init -c genome_config.yaml`, then you add assets using either `refgenie pull` or `refgenie build`. You can also add your own custom assets with `refgenie add`, which is explained in [using custom assets](custom_assets.md).  Refgenie will use the config file to remember what assets are available and where they are.

But here's how the config file works, in case for some reason you do need to edit some things by hand. Here's an example file to get you started: 

```yaml
genome_folder: /path/to/active/genomes
genome_server: http://refgenomes.databio.org
genome_archive: /path/to/archived/genomes
genomes:
    rCRSd:
        assets:
          fasta:
            tags:
              default:
                seek_keys:
                  fasta: rCRSd.fa
                  fai: rCRSd.fa.fai
                  chrom_sizes: rCRSd.chrom.sizes
                asset_parents: []
                asset_path: fasta
                asset_digest: a3c46f201a3ce7831d85cf4a125aa334
                asset_children: ['bowtie2_index:default']
            default_tag: default
          bowtie2_index:
            asset_description: Genome index for bowtie, produced with bowtie-build
            tags:
              default:
                asset_path: bowtie2_index
                seek_keys:
                  bowtie2_index: .
                asset_digest: 0f9217d44264ae2188fcf8e469d2bb38
                asset_parents: ['fasta:default']
            default_tag: default
    hg38:
        assets:
          gencode_gtf:
            asset_description: GTF annotation asset which provides access to all annotated transcripts which make up an Ensembl gene set.
            tags:
              default:
                asset_path: gencode_gtf
                seek_keys:
                  gencode_gtf: hg38.gtf.gz
                asset_digest: 4cd4eac99cdfdeb8ff72d8f8a4a16f9f
                asset_parents: []
            default_tag: default
```

## Details of config attributes

- **genome_folder**: Path to parent folder refgenie-managed assets.
- **genome_server**: URL to a refgenieserver instance (currently only 1 URL allowed).
- **genome_archive**: (optional; used by refgenieserver) Path to folder where asset archives will be stored.
- **genomes**: A list of genomes, each genome has a list of assets. Any relative paths in the asset `path` attributes are considered relative to the genome folder in the config file (or the file itself if not folder path is specified), with the genome name as an intervening path component, e.g. `folder/mm10/indexed_bowtie2`.
- **tags**: A collection of tags defined for the asset
- **default_tag**: A pointer to the tag that is currently defined as the default one
- **asset_parents**: A list of assets that were used to build the asset in question
- **asset_children**: A list of assets that required the asset in question for building
- **seek_keys**: A mapping of names and paths of the specific files within an asset
- **asset_path**: A path to the asset folder, relative to the genome config file
- **asset_digest**: A digest of the asset directory (more precisely, of the file contents within one) used to address the asset provenance issues when the assets are pulled or built.

Note that for a fully operational config just `genome_folder`, `genome_server`, `genomes`, `assets`, `tags` and `seek_keys` keys are required.

For genomes that are managed by `refgenie` (that is, they were built or pulled with `refgenie`), these asset attributes will be automatically populated. You can edit them and refgenie will respect your edits (unless you re-build or re-pull the asset, which will overwrite those fields). You can also add your own assets and `refgenie` won't touch them. For more info, see [using custom assets](custom_assets.md).

## Genome config versions

### v0.2
Up to version `0.4.4`, refgenie used a config file version that lacked the `assets` level in the hierarchy (so, assets were listed directly under the genome). Starting with version `0.5.0`, we moved the assets down a layer to accommodate other genome-level attributes we intend to use in the future (like a description, checksums, other provenance information). Earlier refgenie config files will need to be updated. 

### v0.3
Upt to version `0.6.0`, refgenie used the config v0.2. Currently, it uses v0.3, where we introduced: seek keys, tags, asset digests, default tag pointer, asset description
  
 
