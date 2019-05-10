# Download pre-indexed reference genomes

Most users won't need to become familiar with refgenie because you can just download the pre-indexed archives and use those. These archives are built for common genomes. These are `tar gzipped` files, so **you will need to unarchive them after downloading**. The complete collection is listed at [http://big.databio.org/refgenomes/](http://big.databio.org/refgenomes/):

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