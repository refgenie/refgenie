# Download pre-indexed reference genomes

## Latest builds

The latest genome assets are now available for modular download, either by using the `refgenie` command-line interface, or by browsing on the web. You can browse and download files manually at [refgenomes.databio.org](http://refgenomes.databio.org), but `refgenie` will download and manage these for you automatically. You do this by simply running `refgenie` from the command line.

List all available genomes and assets:

Access [http://refgenomes.databio.org](http://refgenomes.databio.org/assets) or use the CLI to *list remote assets*:

```
refgenie listr
```

Now, you can download the specific asset of your choice with:

```
refgenie pull -g GENOME -a ASSET
```
Where `GENOME` refers to a genome key (*e.g.* hg38) and `ASSET` refers to one or more specific asset keys (*e.g.* bowtie2_index).


## Older builds

For earlier versions of refgenie, there was no `refgenie pull` command, so users would just download the pre-indexed archives and use those. These archives are built for common genomes. These are `tar gzipped` files, so **you will need to unarchive them after downloading**. The complete collection is listed at [http://big.databio.org/refgenomes/](http://big.databio.org/refgenomes/):

## Mirror 1:

* Human: [hg38.tgz](http://big.databio.org/refgenomes/hg38.tgz), [hg19.tgz](http://big.databio.org/refgenomes/hg19.tgz)
* Rat: [rn6.tgz](http://big.databio.org/refgenomes/rn6.tgz)
* Mouse: [mm10](http://big.databio.org/refgenomes/mm10.tgz)
* Prealignment 'decoy' references: available at [big.databio.org/refgenomes/](http://big.databio.org/refgenomes/) for sequences from [ref_decoy](https://github.com/databio/ref_decoy)

## Mirror 2 (use if mirror 1 is down):

* Human: [hg38.tgz](http://cloud.databio.org/refgenomes/hg38.tgz), [hg19.tgz](http://cloud.databio.org/refgenomes/hg19.tgz)
* Rat: [rn6.tgz](http://cloud.databio.org/refgenomes/rn6.tgz)
* Mouse: mm9 (pending), [mm10](http://cloud.databio.org/refgenomes/mm10.tgz)


## Example to download and unarchive hg38_chr22

```
wget http://big.databio.org/refgenomes/hg38_chr22.tgz
tar -xf hg38_chr22.tgz
```