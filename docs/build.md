# Building genome indexes with refgenie

Once you've [installed refgenie](install.md), you can use `refgenie pull` to [download pre-built assets](download.md) without installing any additional software. If you want to build your own, you'll also need to install the building software for the asset you want to build. You have two choices to get that software, you can either [install building software natively](#install_building_software_natively), or use a [docker image](#docker).

Once you're set up with all the additional software, you simply run `refgenie build`, passing it any necessary input files called for by the asset recipe. Further documentation on building specific assets is forthcoming.

## Install building software natively

Refgenie expects to find in your `PATH` any tools needed for building a desired asset. You'll need to follow the instructions for each of these individually. You could find some basic ideas for how to install these programatically in the [dockerfile](https://github.com/databio/refgenie/blob/dev/containers/Dockerfile_refgenie).

Refgenie knows how to build indexes for bowtie2, hisat2, bismark, and other common tools. You can find the complete list in the [config file](https://github.com/databio/refgenie/blob/dev/refgenie/refgenie.yaml). These are all optional; you only have to build indexes for ones you intend to use. You can also add more later. If you don't pass along a configuration file to `refgenie`, it will simply use that one, building those indexes. If you want to choose a subset, copy the config file, edit it as desired, and pass it to `refgenie` like this:
```
wget https://raw.githubusercontent.com/databio/refgenie/master/refgenie/refgenie.yaml
refgenie build -c refgenie.yaml
```

### Adding GTFs

Refgenie also allows you to add information in the form of a GTF file, which provides gene annotation.


## Docker

If you don't want to install all those indexers (and I don't blame you), then you may be interested in my docker image on DockerHub (nsheff/refgenie) that has all of these packages pre-installed, so you can run the complete indexer without worrying about paths and packages. Just clone this repo and run it with the `-d` flag. For example:

```
refgenie build --input rn6.fa -d
```

### Building the container

You can build the docker container yourself like this:

```
git clone https://github.com/databio/refgenie.git
cd refgenie/containers
make refgenie
```

### Pulling the container

```
docker pull nsheff/refgenie
```

# Refgenie Recipes

Here are a few easy scripts you can use to re-index some of your favorite genomes

## hg19

```console
INPUT=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
GTF=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
refgenie build -i ${INPUT} -a ${GTF} -n hg19
```

## hg38
(use the NCBI's official version for sequence alignments without _alt sequences:)
Old link: INPUT=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

This README describes the sequences: 

ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/README_analysis_sets.txt

```console
BUILDER=${CODEBASE}refgenie/src/refgenie.py
INPUT=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
GTF=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_23/gencode.v23.primary_assembly.annotation.gtf.gz
refgenie build -i ${INPUT} -a ${GTF} -n hg38
```

## mm10

```console
BUILDER=${CODEBASE}refgenie/src/refgenie.py
INPUT=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.5_GRCm38.p3/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna.gz
GTF=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.primary_assembly.annotation.gtf.gz
refgenie build -i ${INPUT} -a ${GTF} -n mm10
```
