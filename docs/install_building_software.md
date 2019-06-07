# Installing Genome indexers

 Refgenie expects to find in your PATH any tools needed for building the asset you want to build. You'll need to follow the instructions for each of these individually. You could find some basic ideas for how to install these programatically in the [dockerfile](https://github.com/databio/refgenie/blob/dev/containers/Dockerfile_refgenie).

## List of indexers

You can find a list of indexers you may choose from in the [refgenie config file](https://github.com/databio/refgenie/blob/dev/refgenie/refgenie.yaml).

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

