# Building assets for custom genomes

Once you've [installed refgenie](install.md), you can use `refgenie pull` to [download pre-built assets](download.md) without installing any additional software. However, you may need to use the `build` function for genomes or assets that are not available on the server. You can build assets for any genome for which you can provide the required inputs.

Building assets is a bit more complicated than pulling them. If you want to build assets, you'll need to get the software required by the asset you want to build. You have two choices to get that software: you can either install it natively, or use a docker image (details further down this page). Once you're set up, you simply run `refgenie build`, passing it any necessary input arguments called for by the asset recipe.

## What assets can refgenie build?

At the moment the building functionality is under rapid development and may change in the future. While `refgenie` is totally flexible with respect to genome, it is more restricted in terms of what assets it can build. We are planning to allow users to specify their own recipes for arbitrary assets, but at the moment, `refgenie` can only build a handful of assets for which we have already created building recipes. Refgenie comes with built-in recipes to build indexes for common tools like bowtie2, hisat2, bismark, salmon, bwa, and a few others. If you type `refgenie list`, you'll get a list of all the assets you can build with refgenie. If you want to add a new asset, you'll have to work with us to provide a script that can build it, and we can incorporate it into refgenie. We expect this will get much easier in the future.

Each asset requires some input. For many of the built-in recipes, this is just a FASTA file. Below, we go through the assets you can build and how to build them.

## Examples for top-level assets you can build

### fasta

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

### refgene_anno
A refgene annotation file is used to build several other derived assets.

Some examples:

- hg19: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
- hg38: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
- mm10: http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz
- rn6: http://hgdownload.cse.ucsc.edu/goldenPath/rn6/database/refGene.txt.gz

```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
refgenie build -g hg38 -a refgene_anno --refgene refGene.txt.gz
```

### gencode_gtf

The gencode_gtf asset just copies over a GTF annotation file provided by gencode.

Some examples are:

- hg19: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
- hg38: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gtf.gz
- mm10: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22/gencode.vM22.annotation.gtf.gz

Build the asset like:
```
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
refgenie build -g hg19 -a gencode_gtf --gtf gencode.v19.annotation.gtf.gz
```

### ensembl_gtf
```
wget ftp://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf.gz
refgenie build -g hg38 -a ensembl_gtf --gtf Homo_sapiens.GRCh38.97.gtf.gz
```
Some examples are:

- hg38: ftp://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf.gz
- hg19: ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
- mm10: ftp://ftp.ensembl.org/pub/release-97/gtf/mus_musculus/Mus_musculus.GRCm38.97.gtf.gz
- rn6: ftp://ftp.ensembl.org/pub/release-97/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.97.gtf.gz

### ensembl_rb

This is the ensembl regulatory build. It requires an input `gff` file.

Some examples are:

- hg38: ftp://ftp.ensembl.org/pub/release-96/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190122.gff.gz
- mm10: ftp://ftp.ensembl.org/pub/release-97/regulation/mus_musculus/mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20180516.gff.gz

```
wget ftp://ftp.ensembl.org/pub/release-96/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190122.gff.gz
refgenie build -g hg38 -a ensembl_rb --gff homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190122.gff.gz
```

## Examples for derived assets you can build

### bowtie2 index

The bowtie2_index asset doesn't require any input, but does require that you've already built the `fasta` asset. So, first build the fasta asset for your genome of interest, and then you just build the `bowtie2_index` asset with no other requirements:

```
refgenie build -g test -a bowtie2_index -d 
```

### bismark indexes

The bismark index assets doesn't require any input, but does require that you've already built the `fasta` asset.

```
refgenie build -g test -a bismark_bt2_index -d -R
```

### ensembl_gtf

The ensembl_gtf asset is a copy of the ENSEMBL annotation file. You could build it like this:

```
wget ftp://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf.gz
refgenie build -g hg38 -a ensembl_gtf --gtf Homo_sapiens.GRCh38.97.gtf.gz
```

## Install building software natively

Refgenie expects to find in your `PATH` any tools needed for building a desired asset. You'll need to follow the instructions for each of these individually. You could find some basic ideas for how to install these programatically in the [dockerfile](https://github.com/databio/refgenie/blob/dev/containers/Dockerfile_refgenie). At the moment, the build system is not very flexible, and we don't have great documentation for what is required if you want to use this native approach. In our next major update, we're planning to revamp this system to provide a much more robust build system.

## Building assets with docker

If you don't want to install all the software needed to build all these assets (and I don't blame you), then you can just use docker. Each of our recipes knows about a docker image that has everything it needs. If you have docker installed, you should be able to simply run `refgenie build` with the `-d` flag. For example:

```
refgenie build -d --genome ... --asset ...
```

This tells refgenie to execute the building in a docker container requested by the particular asset recipe you specify. Docker will automatically pull the image it needs when you call this. If you like, you can build the docker container yourself like this:

```
git clone https://github.com/databio/refgenie.git
cd refgenie/containers
make refgenie
```

or pull it directly from [dockerhub](https://hub.docker.com/r/databio/refgenie) like this:

```
docker pull databio/refgenie
```