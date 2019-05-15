# Installing refgenie




## Installing

If you just want to use pre-built refgenie assemblies, just head over to the [download page](download.md); you don't even need to install refgenie. If you want to index your own genomes, then you'll need to install refgenie plus your genome indexers of choice. Install refgenie from [GitHub releases](https://github.com/databio/refgenie/releases) or from PyPI with `pip`:


```console
pip install --user refgenie
```

Update with:

```console
pip install --user --upgrade refgenie
```

See if your install worked by invoking `refgenie` from the command line:

```
refgenie -h
```

If the `refgenie` executable in not automatically in your `$PATH`, add the following line to your `.bashrc` or `.profile` (or `.bash_profile` on MACOS):

```console
export PATH=~/.local/bin:$PATH
```

After that, you'll need to [install the genome indexers](install.md)...

## Genome indexers

Once you've installed refgenie (with `pip install refgenie`), you'll also need the indexers installed. Refgenie expects to find in your PATH any indexer tools you need for the various aligners. You'll need to follow the instructions for each of these individually. You could find some basic ideas for how to install these programatically in the [dockerfile](https://github.com/databio/refgenie/blob/dev/containers/Dockerfile_refgenie).

## List of indexers

You can find a list of indexers you may choose from in the [refgenie config file](https://github.com/databio/refgenie/blob/dev/refgenie/refgenie.yaml).

## Docker

If you don't want to install all those indexers (and I don't blame you), then you may be interested in my docker image on DockerHub (nsheff/refgenie) that has all of these packages pre-installed, so you can run the complete indexer without worrying about paths and packages. Just clone this repo and run it with the `-d` flag. For example:

```
~/code/refgenie/src/refgenie.py --input rn6.fa --outfolder $HOME -d
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

