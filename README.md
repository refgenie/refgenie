# Refgenie: Reference Genome Indexer

Refgenie is a python script that creates a standardized folder structure for reference genome files and indexes. Minimal input is just a fasta file. It also comes with a [docker container (nsheff/refgenie)](https://hub.docker.com/r/nsheff/refgenie/) so you don't have to install any software, if you like.

Refgenie currently builds indexes for tools like bowtie2, hisat2, bismark (for DNA methylation), etc. You can find the complete list in the [config file](src/refgenie.yaml). These are all optional; you only have to build indexes for ones you intend to use. You can also add more later.

## Installing

* Install [Pypiper](http://databio.org/pypiper/) (`pip install --user --upgrade https://github.com/epigen/pypiper/zipball/master`)
* Clone this repo (e.g. `git clone git@github.com:databio/refgenie.git`)
* Install software for indexes to build; put them in your path (default) or specify paths in your [refgenie config file](src/refgenie.yaml). Or, you can use the [Docker version](#docker) and then you don't have to install anything but Docker!

## Running

Run refgenie with: `src/refgenie.py -i INPUT_FILE.fa`. (INPUT_FILE is a fasta file of your reference genome, and can be either a local file or a URL)

#### Optional

* Set an environment shell variable called `GENOMES` to point to where you want your references saved.
* Choose which indexes you want to include by toggling them in the [config file](src/refgenie.yaml).

To build a standard reference for a popular genome, follow one of the [recipes](recipes.md).

## Docker

## Why it's useful

NGS pipelines require reference genomes for alignments and other computation. Usually, pipeline authors have a unique way to organize reference genomes, and they point pipelines to this structure. This means you have to build a new reference genome folder organization for each pipeline, or somehow specify your location to each pipeline. Refgenie standardizes the reference genome index folder structure, so you can build your pipelines around that standardized reference genome format. This makes it easy to switch to a new reference genome, or to share pipelines and reference data among collaborators. Pipelines then need only a single variable: the genome name, and you can even use an environment variable to make all pipelines work seamlessly.

## Other projects

The goal is something like the [iGenomes project](http://support.illumina.com/sequencing/sequencing_software/igenome.html), except they just produce finalized files for you to download, but provide no software to produce your own reference for your own genomes. So, you can't use that to make a standardized format for your internal spike-in genomes or other species they don't provide. Refgenie is a script, so you produce the standard for whatever genome you want.

## Contributing

Pull requests welcome! Add an indexer if you like.
