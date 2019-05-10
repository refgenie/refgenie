# Refgenie Recipes

Here are a few easy scripts you can use to re-index some of your favorite genomes

## hg19

```console
BUILDER=${CODEBASE}refgenie/src/refgenie.py
INPUT=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
GTF=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
${BUILDER} -i ${INPUT} -a ${GTF} -n hg19
```

## hg38
(use the NCBI's official version for sequence alignments without _alt sequences:)
Old link: INPUT=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

This README describes the sequences: 

ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/README_analysis_sets.txt

```console
BUILDER=${CODEBASE}refgenie/src/refgenie.py
INPUT=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
GTF=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_23/gencode.v23.primary_assembly.annotation.gtf.gz
${BUILDER} -i ${INPUT} -a ${GTF} -n hg38
```

## mm10

```console
BUILDER=${CODEBASE}refgenie/src/refgenie.py
INPUT=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.5_GRCm38.p3/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001635.5_GRCm38.p3_no_alt_analysis_set.fna.gz
GTF=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M12/gencode.vM12.primary_assembly.annotation.gtf.gz
${BUILDER} -i ${INPUT} -a ${GTF} -n mm10
```
