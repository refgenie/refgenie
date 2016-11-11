# Reference Genome Indexer (RefGenie) Recipes

## hg19

```
BUILDER=${CODEBASE}refgenie/src/refgenie.py
INPUT=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
GTF=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
${BUILDER} -i ${INPUT} -a ${GTF} -n hg19
```
