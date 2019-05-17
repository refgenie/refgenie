# Building genome indexes with Refgenie

Indexing your own reference genome is as easy as 1-2-3:

1. Install refgenie
2. Choose one: either [install genome indexers](/install) or the [docker image](#docker).
3. Run refgenie with: `refgenie -i INPUT_FILE.fa`. (INPUT_FILE is a fasta file of your reference genome, and can be either a local file or a URL)

## Customizing indexes


Refgenie currently builds indexes for bowtie2, hisat2, bismark (for DNA methylation), and others. You can find the complete list in the [config file](https://github.com/databio/refgenie/blob/dev/refgenie/refgenie.yaml). These are all optional; you only have to build indexes for ones you intend to use. You can also add more later. If you don't pass along a configuration file to `refgenie`, it will simply use this one, building these indexes. If you want to toggle some of them, you may choose which indexes you want to include by toggling them. Just duplicate and edit the config file and pass it to refenie like this:

```
wget https://raw.githubusercontent.com/databio/refgenie/master/refgenie/refgenie.yaml
refgenie -c refgenie.yaml
```

## Adding GTFs

Refgenie also allows you to add information in the form of a GTF file, which provides gene annotation.


#### Optional

* Set an environment shell variable called `GENOMES` to point to where you want your references saved.



