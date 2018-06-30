# Refgenie: Reference Genome Indexer

Refgenie creates a standardized folder structure for reference genome files and indexes. You can download pre-built genomes or use the script to build your own for any genome you like. Minimal input is just a fasta file. It also comes with an optional [docker container (nsheff/refgenie)](https://hub.docker.com/r/nsheff/refgenie/) so you don't have to install any software.

## Rationale

NGS pipelines require reference genomes for alignments and other computation. Usually, each pipeline has a unique way to organize reference genomes, and they may require the assembly to be organized in this structure. The problem with that is you must build a new reference genome folder organization for each pipeline, or at the very least specify your organization to each pipeline specifically. This makes it hard to share pipelines across environments and people, as they have different methods for organizing and passing index files to the pipelines.

If everyone built pipelines to follow a standard structure for reference genomes, then pipelines could just take a string describing that genome (e.g. "hg38") and would be able to know how to find the indexes it needs for that genome.

Refgenie standardizes the reference genome index folder structure, so you can build your pipelines around that standardized reference genome format. This makes it easy to switch to a new reference genome, or to share pipelines and reference data among collaborators. Pipelines then need only a single variable: the genome name, and you can even use an environment variable to make all pipelines work seamlessly. The goal is something like the [iGenomes project](http://support.illumina.com/sequencing/sequencing_software/igenome.html), except they just produce finalized files for you to download, but provide no software to produce your own reference for your own genomes. So, you can't use that to make a standardized format for your internal spike-in genomes or other species they don't provide.

You can download pre-computed tarballs of refgenie assemblies if you like, but more importantly, `refgenie` is a script, so you can produce the standard for whatever genome you want.

## Download pre-built indexed reference genomes

These are built indexes for common genomes. These pre-built downloads are `tar gzipped` files, so **you will need to unarchive them after downloading**. The complete collection is listed at [http://big.databio.org/refgenomes/](http://big.databio.org/refgenomes/):

Mirror 1:
* Human: [hg38.tgz](http://big.databio.org/refgenomes/hg38.tgz), [hg19.tgz](http://big.databio.org/refgenomes/hg19.tgz)
* Rat: [rn6.tgz](http://big.databio.org/refgenomes/rn6.tgz)
* Mouse: [mm10](http://big.databio.org/refgenomes/mm10.tgz)
* Prealignment 'decoy' references: available at [big.databio.org/refgenomes/](http://big.databio.org/refgenomes/) for sequences from [ref_decoy](https://github.com/databio/ref_decoy))

Mirror 2 (use if mirror 1 is down):
* Human: [hg38.tgz](http://cloud.databio.org/refgenomes/hg38.tgz), [hg19.tgz](http://cloud.databio.org/refgenomes/hg19.tgz)
* Rat: [rn6.tgz](http://cloud.databio.org/refgenomes/rn6.tgz)
* Mouse: mm9 (pending), [mm10](http://cloud.databio.org/refgenomes/mm10.tgz)


## Index list

Refgenie currently builds indexes for tools like bowtie2, hisat2, bismark (for DNA methylation), etc. You can find the complete list in the [config file](src/refgenie.yaml). These are all optional; you only have to build indexes for ones you intend to use. You can also add more later.

## Indexing your own reference genome

* Install [Pypiper](http://databio.org/pypiper/) (`pip install --user --upgrade https://github.com/epigen/pypiper/zipball/master`) (Refgenie requires version >= 0.5)
* Clone this repo (e.g. `git clone git@github.com:databio/refgenie.git`)
* Install software for indexes to build; put them in your path (default) or specify paths in your [refgenie config file](src/refgenie.yaml). Or, you can use the [Docker version](#docker) and then you don't have to install anything but [pypiper](http://databio.org/pypiper/) and [docker](http://www.docker.com).

Run refgenie with: `src/refgenie.py -i INPUT_FILE.fa`. (INPUT_FILE is a fasta file of your reference genome, and can be either a local file or a URL)

#### Optional

* Set an environment shell variable called `GENOMES` to point to where you want your references saved.
* Choose which indexes you want to include by toggling them in the [config file](src/refgenie.yaml).

To build a standard reference for a popular genome, follow one of the [recipes](recipes.md).

## Docker

I have produced a docker image on DockerHub (nsheff/refgenie) that has all of these packages pre-installed, so you can run the complete indexer without worrying about paths and packages. Just clone this repo and run it with the `-d` flag. For example:

```
~/code/refgenie/src/refgenie.py -d --input rn6.fa --outfolder $HOME
```

## Contributing

Pull requests welcome! Add an indexer if you like.
