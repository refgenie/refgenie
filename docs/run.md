
To build a standard reference for a popular genome, follow one of the [recipes](recipes.md).

http://big.databio.org/example_data/rCRS.fa.gz
```
refgenie -i http://big.databio.org/example_data/rCRS.fa.gz -R -d
```

## Indexing your own reference genome

* Install [Pypiper](http://databio.org/pypiper/) (`pip install --user --upgrade https://github.com/epigen/pypiper/zipball/master`) (Refgenie requires version >= 0.5)
* Clone this repo (e.g. `git clone git@github.com:databio/refgenie.git`)
* Install software for indexes to build; put them in your path (default) or specify paths in your [refgenie config file](src/refgenie.yaml). Or, you can use the [Docker version](#docker) and then you don't have to install anything but [pypiper](http://databio.org/pypiper/) and [docker](http://www.docker.com).

Run refgenie with: `src/refgenie.py -i INPUT_FILE.fa`. (INPUT_FILE is a fasta file of your reference genome, and can be either a local file or a URL)