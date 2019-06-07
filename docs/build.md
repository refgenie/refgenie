# Building genome indexes with refgenie

Once you've [installed refgenie](/install), you can use `refgenie pull` to [download pre-built assets](/download) without installing any additional software. If you want to build your own, you'll also need to install the building software for the asset you want to build. You have two choices to get that software, you can either [install building software natively](#install_building_software_natively), or use a [docker image](#docker).

Once you're set up with all the additional software, you simply run `refgenie build`, passing it any necessary input files called for by the asset recipe. Further documentation on building specific assets is forthcoming.

## Install building software natively

Refgenie expects to find in your PATH any tools needed for building the asset you want to build. You'll need to follow the instructions for each of these individually. You could find some basic ideas for how to install these programatically in the [dockerfile](https://github.com/databio/refgenie/blob/dev/containers/Dockerfile_refgenie).

Refgenie knows how to build indexes for bowtie2, hisat2, bismark, and others. You can find the complete list in the [config file](https://github.com/databio/refgenie/blob/dev/refgenie/refgenie.yaml). These are all optional; you only have to build indexes for ones you intend to use. You can also add more later. If you don't pass along a configuration file to `refgenie`, it will simply use this one, building these indexes. If you want to toggle some of them, you may choose which indexes you want to include by toggling them. Just duplicate and edit the config file and pass it to refenie like this:

```
wget https://raw.githubusercontent.com/databio/refgenie/master/refgenie/refgenie.yaml
refgenie build -c refgenie.yaml
```

### Adding GTFs

Refgenie also allows you to add information in the form of a GTF file, which provides gene annotation.


## Docker

If you don't want to install all those indexers (and I don't blame you), then you may be interested in my docker image on DockerHub (nsheff/refgenie) that has all of these packages pre-installed, so you can run the complete indexer without worrying about paths and packages. Just clone this repo and run it with the `-d` flag. For example:

```
~/code/refgenie/refgenie/refgenie.py --input rn6.fa --outfolder $HOME -d
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

