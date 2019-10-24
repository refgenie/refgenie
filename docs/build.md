# Build assets
## Introduction
Once you've [installed refgenie](install.md), you can use `refgenie pull` to [download pre-built assets](pull.md) without installing any additional software. However, you may need to use the `build` function for genomes or assets that are not available on the server. You can build assets for any genome for which you can provide the required inputs.

Building assets is a bit more complicated than pulling them. If you want to build assets, you'll need to get the software required by the asset you want to build. You have three choices to get that software: you can either [install it natively](build.md#install-building-software-natively), [use a docker image](build.md#building-assets-with-docker), or you can [use our new bulker manifest](http://bulker.databio.org/en/latest/). This will start a pipeline that will create the requested asset and populate the genome config file for you. You can see the [example build output](build_output.md).  

## Build an asset
Once you're set up, you simply run `refgenie build`, passing it any necessary input arguments called for by the asset recipe. Each asset requires some input. For many of the built-in recipes, this is just a FASTA file. To learn what are the required inputs or other asset dependencies, add an `-r` flag to the `refgenie build` command: 

```
$ refgenie build hg38/bowtie2_index -r

'hg38/bowtie2_index' build requirements: 
- assets: fasta
```  
In this case you'll need to build the `fasta` asset for `hg38` genome before building `bowtie2_index`. Notice how 'fasta' appears under `assets` and not under `arguments`. What this means is that to build a bowtie2 index, you do *not* provide a fasta file as an argument, as you might expect. Instead, you *must already have a fasta asset managed by `refgenie`*. One of the advantages of this is that it allows refgenie to keep a record of how you've built your assets, so `refgenie` can remember the link between this bowtie2 asset and the fasta asset, which turns out to be very useful for maintaining provenance of your assets. It also makes it easier to build these kind of derived assets, because you don't actually have to pass any additional arguments to build them.

## What assets can refgenie build?

At the moment the building functionality is under rapid development and may change in the future. While `refgenie` is totally flexible with respect to genome, it is more restricted in terms of what assets it can build. We are planning to allow users to specify their own recipes for arbitrary assets, but at the moment, `refgenie` can only build a handful of assets for which we have already created building recipes. 

If you type `refgenie list`, you'll get a list of all the assets you can build with refgenie (these show up under *recipes*). 

```
$ refgenie list

Local recipes: bismark_bt1_index, bismark_bt2_index, bowtie2_index, bwa_index, dbnsfp, ensembl_gtf, ensembl_rb, epilog_index, fasta, feat_annotation, gencode_gtf, hisat2_index, kallisto_index, refgene_anno, salmon_index, star_index
```

If you want to add a new asset, you'll have to work with us to provide a script that can build it, and we can incorporate it into `refgenie`. If you have assets that cannot be scripted, or you want to add some other custom asset you may [manually add custom assets](custom_assets.md) and still have them managed by `refgenie`. We expect this will get much easier in the future.

Below, we go through the assets you can build and how to build them.

## Top-level assets you can build

### fasta

We recommend for every genome, you first build the `fasta` asset, because it's a starting point for building a lot of other assets. You just have to give a *compressed* fasta file.

Some examples are:

- [hg19 fasta](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz)
- [hg38 fasta](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)
- [mm10 fasta](ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz)
- [rCRS fasta](http://big.databio.org/example_data/rCRS.fa.gz)

```
wget http://big.databio.org/example_data/rCRS.fa.gz
refgenie build rCRS/fasta --fasta rCRS.fa.gz
refgenie seek rCRS/fasta
```

### refgene_anno
A [refGene annotation file](http://varianttools.sourceforge.net/Annotation/RefGene) is used to build several other derived assets include transcription start sites, exons, introns, and premature mRNA sequences.

Some examples:

- [hg19 refGene](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz)
- [hg38 refGene](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz)
- [mm10 refGene](http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz)
- [rn6 refGene](http://hgdownload.cse.ucsc.edu/goldenPath/rn6/database/refGene.txt.gz)

```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
refgenie build hg38/refgene_anno --refgene refGene.txt.gz
```

### gencode_gtf

The [gencode_gtf](ftp://ftp.ebi.ac.uk/pub/databases/gencode/_README.TXT) asset is a GTF file provided by gencode containing all annotated transcripts.

Some examples are:

- [hg19 comprehensive gene annotation](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh37_mapping/gencode.v32lift37.annotation.gtf.gz)
- [hg38 comprehensive gene annotation](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz)
- [mm10 comprehensive gene annotation](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gz)

```
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gz
refgenie build mm10/gencode_gtf --gencode_gtf gencode.vM23.annotation.gtf.gz
```

### ensembl_gtf

The [ensembl gtf](https://useast.ensembl.org/info/genome/genebuild/genome_annotation.html) asset is an Ensembl provided gene annotation file used to build other derived assets.

Some examples are:

- [hg38 ensembl annotations](ftp://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz)
- [hg19 ensembl annotations](ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz)
- [mm10 ensembl annotations](ftp://ftp.ensembl.org/pub/current_gtf/mus_musculus/Mus_musculus.GRCm38.98.gtf.gz)
- [rn6 ensembl annotations](ftp://ftp.ensembl.org/pub/current_gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.98.gtf.gz)

```
wget ftp://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf.gz
refgenie build hg38/ensembl-gtf --ensembl_gtf Homo_sapiens.GRCh38.97.gtf.gz
```

### ensembl_rb

The ensembl rb asset is the [ensembl regulatory build](http://useast.ensembl.org/info/genome/funcgen/regulatory_build.html). It requires a `gff` file as input, and is then utilized to build several derived assets.

Some examples are:

- [hg38 regulatory build](ftp://ftp.ensembl.org/pub/current_regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff.gz)
- [hg19 regulatory build](ftp://ftp.ensembl.org/pub/grch37/current/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20180925.gff.gz)
- [mm10 regulatory build](ftp://ftp.ensembl.org/pub/current_regulation/mus_musculus/mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20180516.gff.gz)

```
wget ftp://ftp.ensembl.org/pub/current_regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff.gz
refgenie build hg38/ensembl_rb --gff homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff.gz
```

### dbnsfp

The [`dbnsfp` asset](http://varianttools.sourceforge.net/Annotation/dbNSFP) is the annotation database for non-synonymous SNPs. It requires the dbNSFP zip file.

```
wget ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFP4.0a.zip
refgenie build test/dbnsfp --dbnsfp dbNSFP4.0a.zip
```

## Derived assets you can build

For many of the following derived assets, you will need the corresponding software to build the asset.  You can either [install software on a case-by-case basis natively](build.md#install-building-software-natively), or you can [build the assets using `docker`](build.md#building-assets-with-docker).

### bowtie2 index

<i class="fas fa-exclamation-circle"></i> requires [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

The `bowtie2_index` asset doesn't require any input, but does require that you've already built the `fasta` asset. So, first [build the `fasta` asset](build.md#fasta) for your genome of interest, and then you just build the `bowtie2_index` asset with no other requirements:

```
refgenie build test/bowtie2_index
```

### bismark indexes

<i class="fas fa-exclamation-circle"></i> requires [bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)


The `bismark_bt1_index` and `bismark_bt2_index`  assets don't require any input, but do require that you've already [built the `fasta` asset](build.md#fasta).

```
refgenie build test/bismark_bt1_index
refgenie build test/bismark_bt2_index
```

### bwa index

<i class="fas fa-exclamation-circle"></i> requires [bwa](http://bio-bwa.sourceforge.net/)

The `bwa_index` asset doesn't require any input, but does require that you've already [built the `fasta` asset](build.md#fasta). Then, you just build the `bwa_index` asset with no other requirements:

```
refgenie build test/bwa_index
```

### hisat2 index

<i class="fas fa-exclamation-circle"></i> requires [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml)

The `hisat2_index` asset doesn't require any input, but does require that you've already [built the `fasta` asset](build.md#fasta). Then, you just build the `hisat2_index` asset with no other requirements:

```
refgenie build test/hisat2_index
```

### epilog index

<i class="fas fa-exclamation-circle"></i> requires [epilog](https://github.com/databio/epilog)

The `epilog_index` asset [requires the call site context as input](https://github.com/databio/epilog), in addition to requiring that you've already [built the `fasta` asset](build.md#fasta). To build an epilog index, you pass the context as an additional argument:

```
refgenie build test/epilog_index --context CG
```

### kallisto index

<i class="fas fa-exclamation-circle"></i> requires [kallisto](https://pachterlab.github.io/kallisto/)

The `kallisto_index` asset doesn't require any input, but does require that you've already [built the `fasta` asset](build.md#fasta). Then, you just build the `kallisto_index` asset with no other requirements:

```
refgenie build test/kallisto_index
```

### salmon index

<i class="fas fa-exclamation-circle"></i> requires [salmon](https://salmon.readthedocs.io/en/latest/salmon.html)

The `salmon_index` asset doesn't require any input, but does require that you've already [built the `fasta` asset](build.md#fasta). Then, you just build the `salmon_index` asset with no other requirements:

```
refgenie build test/salmon_index
```

### star index

<i class="fas fa-exclamation-circle"></i> requires [star](https://github.com/alexdobin/STAR)

The `star_index` asset doesn't require any input, but does require that you've already [built the `fasta` asset](build.md#fasta). Then, you just build the `star_index` asset with no other requirements:

```
refgenie build test/star_index
```

### feature annotation

The `feat_annotation` asset is a combination of genomic feature annotations from [`ensembl_gtf`](build.md#ensembl_gtf) and [`ensembl_rb`](build.md#ensembl_rb) assets.  It requires both assets, but does not require any input as long as both assets have already been built. The resulting asset includes the following genomic feature annotations: enhancers, promoters, promoter flanking regions, 5' UTR, 3' UTR, exons, and introns.

```
refgenie build test/feat_annotation
```

## Install building software natively

`Refgenie` expects to find in your `PATH` any tools needed for building a desired asset. You'll need to follow the instructions for each of these individually. You could find some basic ideas for how to install these programatically in the [dockerfile](https://github.com/databio/refgenie/blob/dev/containers/Dockerfile_refgenie). At the moment, the build system is not very flexible, and we don't have great documentation for what is required if you want to use this native approach. In our next major update, we're planning to revamp this system to provide a much more robust build system.

## Building assets with docker

If you don't want to install all the software needed to build all these assets (and I don't blame you), then you can just use `docker`. Each of our recipes knows about a `docker image` that has everything it needs. If you have `docker` installed, you should be able to simply run `refgenie build` with the `-d` flag. For example:

```
refgenie build -d genome/asset ...
```

This tells `refgenie` to execute the building in a `docker` container requested by the particular asset recipe you specify. `Docker` will automatically pull the image it needs when you call this. If you like, you can build the `docker` container yourself like this:

```
git clone https://github.com/databio/refgenie.git
cd refgenie/containers
make refgenie
```

or pull it directly from [dockerhub](https://hub.docker.com/r/databio/refgenie) like this:

```
docker pull databio/refgenie
```

## Versioning the assets

`refgenie` supports tags to facilitate management of multiple "versions" of the same asset. Simply add a `:your_tag_name` appendix to the asset registry path in the `refgenie build` command and the created asset will be tagged:

```
refgenie build hg38/bowtie2_index:my_tag
```

You can also learn more about [tagging refgenie assets](tag.md).