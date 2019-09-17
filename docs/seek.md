# Retrieve paths to assets

Once you've assembled a few assets, either by downloading or by building them, you'll be able to use `refgenie seek` to retrieve the paths. It's quite simple, really -- say you've built the `bowtie2_index` and `fasta` assets for `hg38`. If you type:

```console
refgenie seek -c CONFIG.yaml hg38/bowtie2_index
```

You'll get back the absolute path on your system to the `bowtie2_index` asset, something like:

```console
/path/to/genomes_folder/hg38/bowtie2_index
```

Because you have also built the `fasta` asset, you'll have available a few more asset keys to retrieve: `fai`, `fasta`, and `chrom_sizes`. Why is this useful? Instead of writing your shell script to take a path to each the `chrom_sizes` file and the `bowtie2_index` path, you just take the genome name, and use `refgenie seek` to find the correct path. Not only have you simplified the interface to your analysis, but both the script itself *and your call of it* are now portable.


```console
GENOME="hg38"

refgenie seek -g $GENOME -a bowtie2_index
refgenie seek -g $GENOME -a fasta.chrom_sizes
refgenie seek -g $GENOME -a fasta
```



