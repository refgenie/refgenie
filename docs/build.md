# Building genome indexes with refgenie

Once you've [installed refgenie](install.md), you can use `refgenie pull` to [download pre-built assets](download.md) without installing any additional software. However, you may need to use the `build` function for genomes or assets that are not available on the server. 

If you want to build assets, you'll need to install the building software for the asset you want to build. You have two choices to get that software: you can either install building software natively, or use a docker image. Once you're set up, you simply run `refgenie build`, passing it any necessary input arguments called for by the asset recipe.

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

## Building assets for custom genomes

You can build assets for any genome for which you can provide the required inputs. For a typical genome index, this would be just a FASTA file, or perhaps another annotation file.

## Building custom assets

At the moment the building functionality is under rapid development and may change in the future. While `refgenie` is totally flexible with respect to genome, it is more restricted in terms of what assets it can build. We are planning to allow users to specify their own recipes for arbitrary assets, but at the moment, `refgenie` can only build a handful of assets for which we have already created building recipes. If you want to add a new asset, you'll have to work with us to provide a script that can build it, and we can incorporate it into refgenie. We expect this will get much easier in the future.


# Examples for assets you can build

## fasta

We recommend for every genome, you first build the `fasta` asset, because it's a starting point for building a lot of other assets. You just have to give a compressed fasta file.

Some examples are:

- hg19: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
- hg38: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
- mm10: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.5_GRCm38.p3/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna.gz
- This [README](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/README_analysis_sets.txt) describes the sequences.


```
export REFGENIE="test.yaml"
refgenie build -g test -a fasta --fasta rCRS.fa.gz
refgenie seek -g test -a fasta
```

## bowtie2 index

Once you have the fasta asset built, you can pass it automatically like this:

```
refgenie build -g test -a bowtie2_index --fasta $(refgenie seek -g test -a fasta) -d -R
```

## bismark indexes

The bismark asset doesn't require any input, but does require that you've already built the `fasta` asset.

```
refgenie build -g test -a bismark_bt2_index -d -R
```

## gtf_anno

The get_anno asset just copies over a GTF annotation file, such as one provided by gencode.

- hg19: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
- hg38: GTF=ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.primary_assembly.annotation.gtf.gz
- mm10: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.primary_assembly.annotation.gtf.gz
Build the asset like:
```
refgenie build -g hg19 -a gtf_anno --gtf ${GTF}
```

