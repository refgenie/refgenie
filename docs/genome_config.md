# Refgenie genome configuration file

Refgenie will read and write a genome configuration file in yaml format. In general, you shouldn't need to mess with the config file. You create one with `refgenie init -c genome_config.yaml`, then you add assets using either `refgenie pull` or `refgenie build`. You can also add your own custom assets with `refgenie add`, which is explained in [using external assets](external_assets.md).  Refgenie will use the config file to remember what assets are available and where they are.

But here's how the config file works, in case for some reason you do need to edit some things by hand. Here's an example file to get you started: 

```yaml
genome_folder: /path/to/active/genomes
genome_server: http://refgenomes.databio.org
genome_archive: /path/to/archived/genomes
genomes:
  hg38:
    bowtie2_index:
      path: indexed_bowtie2
      asset_description: ...
    hisat2_index: 
      path: indexed_hisat2
      asset_description: ...
  mm10:
      bowtie2_index:
        path: indexed_bowtie2
  rCRS:
    bowtie2_index:
      path: indexed_bowtie2
    indexed_bowtie2:
      path: indexed_bowtie2
    indexed_hisat2:
      path: indexed_hisat2
    indexed_kallisto:
      path: indexed_kallisto

```

## Details of config attributes

- **genome_folder**: Path to parent folder refgenie-managed assets.
- **genome_server**: URL to a refgenieserver instance (currently only 1 URL allowed).
- **genome_archive**: (optional; used by refgenieserver) Path to folder where asset archives will be stored.
 - **genomes**: A list of genomes, each genome has a list of assets. Any relative paths in the asset `path` attributes are considered relative to the genome folder in the config file (or the file itself if not folder path is specified), with the genome name as an intervening path component, e.g. `folder/mm10/indexed_bowtie2`.

So, for example:

```yaml
genomes:
  hg38:                <---- genome key
    bowtie2:           <---- asset key
      path: blahblah   <---- relative to genome folder
      asset_description: freeform text (currently manual)
      asset_checksum: 12345
      asset_size: 15MB
      archive_size: 12MB
```

For genomes that are managed by refgenie (that is, they were built or pulled with `refgenie`), these asset attributes will be automatically populated. You can edit them and refgenie will respect your edits (unless you re-build or re-pull the asset, which will overwrite those fields). You can also add your own assets and refgenie won't touch them. For more info, see [using external assets](external_assets.md).
