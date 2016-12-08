# Refgenie: Reference Genome Indexer

Refgenie is a python script that creates a standardized folder structure for reference genome files and indexes. Minimal input is just a fasta file. 

## Why it's useful

NGS pipelines require reference genomes for alignments and other computation. Usually, pipeline authors have a unique way to organize reference genomes, and they point pipelines to this structure. This means you have to build a new reference genome folder organization for each pipeline, or somehow specify your location to each pipeline. Refgenie standardizes the reference genome index folder structure, so you can build your pipelines around that standardized reference genome format. This makes it easy to switch to a new reference genome, or to share pipelines and reference data among collaborators. Pipelines then need only a single variable: the genome name, and you can even use an environment variable to make all pipelines work seamlessly.

Refgenie currently builds indexes for:
* bowtie
* bowtie2
* STAR
* bismark (bt1 and bt2)
* kallisto

## Prerequisites

You need: 
* [Pypiper](http://databio.org/pypiper/) (`pip install --user https://github.com/epigen/pypiper/zipball/master`)
* Indexers for any indexes you want to build; put them in your path or specify them in the [config file](src/refgenie.yaml)

## How to use

1. Clone this repo
2. Run refgenie with: `src/refgenie.py -i INPUT_FILE.fa`. (INPUT_FILE is a fasta file of your reference genome, and can be either a local file or a URL)

#### Optional
* Set an environment shell variable called `GENOMES` to point to where you want your references saved.
* Choose which indexes you want to include by toggling them in the [config file](src/refgenie.yaml).

To build a standard reference for a popular genome, follow one of the [recipes](recipes.md).

## Other projects

The goal is something like the [iGenomes project](http://support.illumina.com/sequencing/sequencing_software/igenome.html), except they just produce finalized files for you to download, but provide no software to produce your own reference for your own genomes. So, you can't use that to make a standardized format for your internal spike-in genomes or other species they don't provide. Refgenie is a script, so you produce the standard for whatever genome you want.

## Contributing

Pull requests welcome!
