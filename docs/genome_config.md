# Refgenie genome configuration file

Refgenie will read and write a genome configuration file in yaml format.

In the future, you shouldn't need to really worry about it... you use `refgenie init` to create a file, then `refgenie pull` and `refgenie build` will populate it as necessary.

But for now we need to document it and you might have to (or want to) edit some things by hand. Here's an example file to get you started: 

```yaml
genome_folder: /path/to/active/genomes
genome_server: http://refgenomes.databio.org
genome_archive: /path/to/archived/genomes
genomes:
  hg38:
    bowtie2: indexed_bowtie2
    hisat2: indexed_hisat2
    indexed_bowtie2: indexed_bowtie2
  mm10:
    indexes:
      bowtie2: indexed_bowtie2
  rCRS:
    bowtie2_index:
      path: indexed_bowtie2
    indexed_bowtie2:
      path: rCRS/indexed_bowtie2
    indexed_hisat2:
      path: rCRS/indexed_hisat2
    indexed_kallisto:
      path: rCRS/indexed_kallisto

```

## Top-level attributes

### genome_folder

### genome_server

### genome_archive

### genomes

Here you have a list of genomes, each genome has a list of assets...

```yaml
genomes:
  genome_name:
  	asset1:
  	  path: path_to_asset
  	  asset_description: ...

```

Any relative paths in the asset `path` attributes are considered relative to the genome config file.

So, for example:

```
genomes:
  hg38:                <--------------- that's the genome key.
    bowtie2:           <--------------- that's the asset key.
      path: blahblah   <--------------- relative to genome config. returned by get_asset('hg38', 'bowtie2'). 
```
