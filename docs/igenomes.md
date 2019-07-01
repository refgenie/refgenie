# Using refgenie with iGenomes

If you're already using iGenomes, it's easy to configure refgenie to use your existing folder structure. [iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html) is project that provides sequences and annotation files for commonly analyzed organisms. Each iGenome is available as a compressed file that contains sequences and annotation files for a single genomic build of an organism. 

Initialize a refgenie config file if you don't have one you want to use for your iGenomes assets:

```
export REFGENIE='genome_config.yaml'
refgenie init -c $REFGENIE
```

And then add individual assets you want refgenie to track with `refgenie add`:

```
refgenie add -g GENOME -a ASSET -p RELATIVE_PATH
```

So, for example, 

```console
wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Caenorhabditis_elegans/UCSC/ce10/Caenorhabditis_elegans_UCSC_ce10.tar.gz
tar -xf Caenorhabditis_elegans_UCSC_ce10.tar.gz
refgenie init -c igenome_config.yaml
refgenie add -c igenome_config.yaml -g ce10 --asset bowtie2_index --path Caenorhabditis_elegans/UCSC/ce10/Sequence/Bowtie2Index
```

Now we can `seek` any of those assets:

```console
refgenie seek -c igenome_config.yaml -g ce10 --asset bowtie2_index
```

This way you can configure refgenie to use your iGenomes assets, so you can wean yourself off of the iGenomes hard structure and transition to the refgenie-managed path system.
