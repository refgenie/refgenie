# Building genome indexes with refgenie

Once you've [installed refgenie](install.md), you can use `refgenie pull` to [download pre-built assets](download.md) without installing any additional software. If you want to build assets instead of downloading them, you'll also need to install the building software for the asset you want to build. You have two choices to get that software, you can either [install building software natively](#install_building_software_natively), or use a [docker image](#docker).

Once you're set up with all the additional software, you simply run `refgenie build`, passing it any necessary input arguments called for by the asset recipe.

## Install building software natively

Refgenie expects to find in your `PATH` any tools needed for building a desired asset. You'll need to follow the instructions for each of these individually. You could find some basic ideas for how to install these programatically in the [dockerfile](https://github.com/databio/refgenie/blob/dev/containers/Dockerfile_refgenie).

Refgenie knows how to build indexes for bowtie2, hisat2, bismark, and other common tools. You can find the complete list in the [config file](https://github.com/databio/refgenie/blob/dev/refgenie/refgenie.yaml).

## Building assets with docker

If you don't want to install all the software needed to build all these assets (and I don't blame you), then you may be interested in our docker image on DockerHub (nsheff/refgenie) that has all of these packages pre-installed, so you can run the complete indexer without worrying about paths and packages. 

You can build the docker container yourself like this:

```
git clone https://github.com/databio/refgenie.git
cd refgenie/containers
make refgenie
```

or pull it directly from dockerhub like this:

```
docker pull nsheff/refgenie
```

Once you have the container, run `refgenie` with the `-d` flag. For example:

```
refgenie build -d --genome ... --asset ...
```

# Refgenie Recipes

Here are a few easy scripts you can use to re-index some of your favorite genomes

## hg19

```console
FASTA=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit

wget ${FASTA}
twoBitToFa hg19.2bit hg19.fa
refgenie build -g hg19 -a fasta --fasta hg19.fa
```

```
GTF=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
refgenie build -g hg19 -a gtf_anno --gtf ${GTF}
```

## hg38
(use the NCBI's official version for sequence alignments without _alt sequences:)

Old link: INPUT=`ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz`

This README describes the sequences: 

[ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/README_analysis_sets.txt](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/README_analysis_sets.txt)

```console
INPUT=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
GTF=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.primary_assembly.annotation.gtf.gz
refgenie build -i ${INPUT} -a ${GTF} -n hg38
```

## mm10

```console
INPUT=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.5_GRCm38.p3/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna.gz
GTF=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.primary_assembly.annotation.gtf.gz
refgenie build -i ${INPUT} -a ${GTF} -n mm10
```
