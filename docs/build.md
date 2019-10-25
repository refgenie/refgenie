# Build assets
## Introduction
Once you've [installed refgenie](install.md), you can use `refgenie pull` to [download pre-built assets](pull.md) without installing any additional software. However, you may need to use the `build` function for genomes or assets that are not available on the server. You can build assets for any genome for which you can provide the required inputs.

Building assets is a bit more complicated than pulling them. You'll need to set up 2 things: 1) any *inputs* required by the asset build recipe; 2) Any *software* required by the recipe. Below, we'll walk you through each of these requirements, but first, how can you tell *what* refgenie can build in the first place?

## What assets can refgenie build?

At the moment the building functionality is under rapid development and may change in the future. While `refgenie` is totally flexible with respect to genome, it is more restricted in terms of what assets it can build. We are planning to allow users to specify their own recipes for arbitrary assets, but at the moment, `refgenie` can only build a handful of assets for which we have already created building recipes. If you type `refgenie list`, you'll get a list of all the assets you can build with refgenie (under *recipes*). You can also browse the [list of available assets](available_assets.md) here. If you need refgenie to manage an asset not in this list, you can either 1) wait for our pending implementation of custom recipes, or 2 [add custom assets](custom_assets.md), which you would build separately and then use refgenie just to manage them.


## 1. Required asset inputs

Each asset requires some inputs, which can be either arguments (external files) or pre-existing assets already managed by refgenie. To view required inputs for an asset, add an `-r` flag to the `refgenie build` command: 

```
$ refgenie build hg38/bowtie2_index -r

'hg38/bowtie2_index' build requirements: 
- assets: fasta
```  
Notice how 'fasta' appears under `assets` and not under `arguments`. This means to build a bowtie2 index, you do *not* provide a fasta file as an *argument*; instead, you *must already have a fasta asset managed by `refgenie`*. One advantage of this is that it allows refgenie to keep a record of how you've built your assets, so `refgenie` can remember the link between this bowtie2 asset and the fasta asset, which turns out to be very useful for maintaining provenance of your assets. It also makes it easier to build derived assets like this, because you don't actually have to pass any additional arguments to build them.

So, you'll need to build the `fasta` asset for `hg38` genome before building `bowtie2_index`, but once you have that, building this asset is as simple as typing:

```
$ refgenie build hg38/bowtie2_index
```

For many of the built-in recipes, a pre-existing FASTA asset is the only requirement. Next, here's an example of an asset that requires an argument, but not a pre-existing asset:

```
$ refgenie build hg38/refgene_anno -r

'hg38/refgene_anno' build requirements: 
- arguments: refgene
```

You'll need to provide this recipe with a `refgene` argument, like this:

```
$ refgenie build hg38/refgene_anno --refgene REFGENE_FILE.txt.gz
```

You can see the [example build output](build_output.md).

## 2. Required asset software

If you want to build assets, you'll need to get the software required by the asset you want to build. You have three choices to get that software: you can either install it natively, use a docker image, or use a bulker manifest.   

### Install build software natively

`Refgenie` expects to find in your `PATH` any tools needed for building a desired asset. You'll need to follow the instructions for each of these individually. You could find some basic ideas for how to install these programatically in the [dockerfile](https://github.com/databio/refgenie/blob/dev/containers/Dockerfile_refgenie). We discourage this approach because it makes the assets dependent on your particular uncontrolled environment, which is not ideal. As a result, we don't have great documentation for what is required if you want to use this native approach. As we develop a custom asset system, we're planning to revamp this to provide more detailed way to see what requirements are for a specific recipe.

### Build assets with docker

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

### Build assets with bulker

For an even more seamless integration of containers with `refgenie`, learn about [`bulker`](http://bulker.databio.org/en/latest/refgenie_tutorial/), our multi-container environment manager. Here, you'd just need to do this:

```
pip install bulker

# Next, configure bulker according to your local compute environment

bulker load databio/refgenie:0.7.0
bulker activate databio/refgenie:0.7.0
refgenie build ...
```

Bulker works on both singularity and docker systems.

## Versioning the assets

`refgenie` supports tags to facilitate management of multiple "versions" of the same asset. Simply add a `:your_tag_name` appendix to the asset registry path in the `refgenie build` command and the created asset will be tagged:

```
refgenie build hg38/bowtie2_index:my_tag
```

You can also learn more about [tagging refgenie assets](tag.md).