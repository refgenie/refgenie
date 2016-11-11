Refgenie: Reference Genome Indexer

Refgenie is a scripted way to create a standardized format for reference genome files and indexes. It's something like the [iGenomes project](http://support.illumina.com/sequencing/sequencing_software/igenome.html), except they just produce finalized files for you to download, but provide no software to produce your own reference for your own genomes. So, you can't use that to make a standardized format for your internal spike-in genomes or other species they don't provide. Refgenie is a script, so you produce the standard for whatever genome you want.

It currently builds indexes for:
* bowtie
* bowtie2
* STAR
* bismark (bt1 and bt2)
* kallisto

## Prerequisites

You need: 
* [Pypiper](http://databio.org/pypiper/)
* Indexers in your path for any indexes you want to build, or, specify them in the [config file](src/refgenie.yaml)

## How to use

1. Clone this repo
2. Set an environment variable `GENOMES` to point to where you want your references saved (optional).
3. Follow one of the [recipes](recipes.md). At a minimum, you just need to provide a fasta file, either local or as a URL.

## Contributing

Pull requests welcome!