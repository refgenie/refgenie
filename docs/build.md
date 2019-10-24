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

At the moment the building functionality is under rapid development and may change in the future. While `refgenie` is totally flexible with respect to genome, it is more restricted in terms of what assets it can build. We are planning to allow users to specify their own recipes for arbitrary assets, but at the moment, `refgenie` can only build a handful of [assets for which we have already created building recipes](available_assets.md). 

If you type `refgenie list`, you'll get a [list of all the assets you can build with refgenie](available_assets.md) (these show up under *recipes*). 

```
$ refgenie list

Local recipes: bismark_bt1_index, bismark_bt2_index, bowtie2_index, bwa_index, dbnsfp, ensembl_gtf, ensembl_rb, epilog_index, fasta, feat_annotation, gencode_gtf, hisat2_index, kallisto_index, refgene_anno, salmon_index, star_index
```

You can also [add custom assets](custom_assets.md) you want `refgenie` to manage.

## Install build software natively

`Refgenie` expects to find in your `PATH` any tools needed for building a desired asset. You'll need to follow the instructions for each of these individually. You could find some basic ideas for how to install these programatically in the [dockerfile](https://github.com/databio/refgenie/blob/dev/containers/Dockerfile_refgenie). At the moment, the build system is not very flexible, and we don't have great documentation for what is required if you want to use this native approach. In our next major update, we're planning to revamp this system to provide a much more robust build system.

## Build assets with docker

If you don't want to install all the software needed to build all these assets (and I don't blame you), then you can just use `docker`. Each of our recipes knows about a `docker image` that has everything it needs. If you have `docker` installed, you should be able to simply run `refgenie build` with the `-d` flag. For example:

```
refgenie build -d genome/asset ...
```

This tells `refgenie` to execute the building in a `docker container` requested by the particular asset recipe you specify. `Docker` will automatically pull the image it needs when you call this. If you like, you can build the `docker container` yourself like this:

```
git clone https://github.com/databio/refgenie.git
cd refgenie/containers
make refgenie
```

or pull it directly from [dockerhub](https://hub.docker.com/r/databio/refgenie) like this:

```
docker pull databio/refgenie
```

For an even more seamless integration of containers with `refgenie`, [learn about `bulker`, our multi-container environment manager](http://bulker.databio.org/en/latest/refgenie_tutorial/).

## Versioning the assets

`refgenie` supports tags to facilitate management of multiple "versions" of the same asset. Simply add a `:your_tag_name` appendix to the asset registry path in the `refgenie build` command and the created asset will be tagged:

```
refgenie build hg38/bowtie2_index:my_tag
```

You can also learn more about [tagging refgenie assets](tag.md).