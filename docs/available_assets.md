<style>
.fas {
  width: 25px;
  margin-right: 5px;
  text-align: center;
  horizontal-align: center;
}
</style>

# Buildable assets

`Refgenie` can build a handful of assets for which we have already created building recipes. `refgenie list` lists all assets refegenie can build:

```
$ refgenie list

Local recipes: bismark_bt1_index, bismark_bt2_index, bowtie2_index, bwa_index, dbnsfp, ensembl_gtf, ensembl_rb, epilog_index, fasta, feat_annotation, gencode_gtf, hisat2_index, kallisto_index, refgene_anno, salmon_index, star_index, suffixerator_index, tallymer_index
```

If you want to add a new asset, you'll have to work with us to provide a script that can build it, and we can incorporate it into `refgenie`. If you have assets that cannot be scripted, or you want to add some other custom asset you may [manually add custom assets](custom_assets.md) and still have them managed by `refgenie`. We expect this will get much easier in the future.

Below, we go through the assets you can build and how to build them.

## Top-level assets you can build

### fasta

<i class="fas fa-file-import"></i> required files: `--files fasta=/path/to/fasta_file` (*e.g.* [example_genome.fa.gz](http://big.databio.org/example_data/rCRS.fa.gz))
<i class="fas fa-sliders-h"></i> required parameters: *none*
<i class="fas fa-exclamation-triangle"></i> required asset: *none*
<i class="fas fa-tools"></i> required software: [samtools](http://www.htslib.org/)

We recommend for every genome, you first build the `fasta` asset, because it's a starting point for building a lot of other assets.

Example fasta files:

- [hg19 fasta](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz)
- [hg38 fasta](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)
- [mm10 fasta](ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz)
- [rCRS fasta](http://big.databio.org/example_data/rCRS.fa.gz)

```
wget http://big.databio.org/example_data/rCRS.fa.gz
refgenie build rCRS/fasta --files fasta=rCRS.fa.gz
refgenie seek rCRS/fasta
```

### blacklist

<i class="fas fa-file-import"></i> required files: `--files blacklist=/path/to/blacklist_file` (*e.g.* [hg38-blacklist.v2.bed.gz](https://github.com/Boyle-Lab/Blacklist/tree/master/lists))
<i class="fas fa-sliders-h"></i> required parameters: *none*
<i class="fas fa-exclamation-triangle"></i> required asset: *none*
<i class="fas fa-tools"></i> required software: *none*

The `blacklist` asset represents regions that should be excluded from sequencing experiments. The ENCODE blacklist represents a comprehensive listing of these regions for several model organisms [^Amemiya2019].

Example blacklist files:

- [hg19 blacklist](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg19-blacklist.v2.bed.gz)
- [hg38 blacklist](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz)
- [mm10 blacklist](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz)
- [dm6 blacklist](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/dm6-blacklist.v2.bed.gz)

```
wget https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg38-blacklist.v2.bed.gz
refgenie build hg38/blacklist --files blacklist=hg38-blacklist.v2.bed.gz
```

### refgene_anno

<i class="fas fa-file-import"></i> required files: `--files refgene=/path/to/refGene_file` (*e.g.* [refGene.txt.gz](http://varianttools.sourceforge.net/Annotation/RefGene))
<i class="fas fa-sliders-h"></i> required parameters: *none*
<i class="fas fa-exclamation-triangle"></i> required asset: *none*
<i class="fas fa-tools"></i> required software: *none*

The `refgene_anno` asset is used to produce derived assets including transcription start sites (TSSs), exons, introns, and premature mRNA sequences.

Example refGene annotation files:

- [hg19 refGene](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz)
- [hg38 refGene](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz)
- [mm10 refGene](http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/refGene.txt.gz)
- [rn6 refGene](http://hgdownload.cse.ucsc.edu/goldenPath/rn6/database/refGene.txt.gz)

```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
refgenie build hg38/refgene_anno --files refgene=refGene.txt.gz
```

### gencode_gtf

<i class="fas fa-file-import"></i> required files: `--files gencode_gtf=/path/to/gencode_file` (*e.g.* [gencode.gtf.gz](ftp://ftp.ebi.ac.uk/pub/databases/gencode/_README.TXT))
<i class="fas fa-sliders-h"></i> required parameters: *none*
<i class="fas fa-exclamation-triangle"></i> required asset: *none*
<i class="fas fa-tools"></i> required software: *none*

The `gencode_gtf` asset contains all annotated transcripts.

Example gencode files:

- [hg19 comprehensive gene annotation](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh37_mapping/gencode.v32lift37.annotation.gtf.gz)
- [hg38 comprehensive gene annotation](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz)
- [mm10 comprehensive gene annotation](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gz)

```
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M23/gencode.vM23.annotation.gtf.gz
refgenie build mm10/gencode_gtf --files gencode_gtf=gencode.vM23.annotation.gtf.gz
```

### ensembl_gtf

<i class="fas fa-file-import"></i> required files: `--files ensembl_gtf=/path/to/ensembl_file` (*e.g.* [ensembl.gtf.gz](https://useast.ensembl.org/info/genome/genebuild/genome_annotation.html))
<i class="fas fa-sliders-h"></i> required parameters: *none*
<i class="fas fa-exclamation-triangle"></i> required asset: *none*
<i class="fas fa-tools"></i> required software: *none*

The `ensembl_gtf` asset is used to build other derived assets including a comprehensive TSS annotation and gene body annotation.

Example Ensembl files:

- [hg38 ensembl annotations](ftp://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz)
- [hg19 ensembl annotations](ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz)
- [mm10 ensembl annotations](ftp://ftp.ensembl.org/pub/current_gtf/mus_musculus/Mus_musculus.GRCm38.98.gtf.gz)
- [rn6 ensembl annotations](ftp://ftp.ensembl.org/pub/current_gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.98.gtf.gz)

```
wget ftp://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf.gz
refgenie build hg38/ensembl-gtf --files ensembl_gtf=Homo_sapiens.GRCh38.97.gtf.gz
```

### ensembl_rb

<i class="fas fa-file-import"></i> required files: `--files gff=/path/to/gff_file` (*e.g.* [regulatory_features.ff.gz](http://useast.ensembl.org/info/genome/funcgen/regulatory_build.html))
<i class="fas fa-sliders-h"></i> required parameters: *none*
<i class="fas fa-exclamation-triangle"></i> required asset: *none*
<i class="fas fa-tools"></i> required software: *none*

The `ensembl_rb` asset is used to produce derived assets including feature annotations.

Example Ensembl files:

- [hg38 regulatory build](ftp://ftp.ensembl.org/pub/current_regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff.gz)
- [hg19 regulatory build](ftp://ftp.ensembl.org/pub/grch37/current/regulation/homo_sapiens/homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20180925.gff.gz)
- [mm10 regulatory build](ftp://ftp.ensembl.org/pub/current_regulation/mus_musculus/mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20180516.gff.gz)

```
wget ftp://ftp.ensembl.org/pub/current_regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff.gz
refgenie build hg38/ensembl_rb --files gff=homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff.gz
```

### dbnsfp

<i class="fas fa-file-import"></i> required files: `--files dbnsfp=/path/to/dbnsfp_file` (*e.g.* [dbNSFP4.0a.zip](http://varianttools.sourceforge.net/Annotation/dbNSFP))
<i class="fas fa-sliders-h"></i> required parameters: *none*
<i class="fas fa-exclamation-triangle"></i> required asset: *none*
<i class="fas fa-tools"></i> required software: *none*

The `dbnsfp` asset is the annotation database for non-synonymous SNPs.

```
wget ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFP4.0a.zip
refgenie build test/dbnsfp --files dbnsfp=dbNSFP4.0a.zip
```

## Derived assets you can build

For many of the following derived assets, you will need the corresponding software to build the asset.  You can either [install software on a case-by-case basis natively](build.md#install-building-software-natively), or you can [build the assets using `docker`](build.md#building-assets-with-docker).

### bowtie2_index

<i class="fas fa-file-import"></i> required files: *none*
<i class="fas fa-sliders-h"></i> required parameters: *none*
<i class="fas fa-exclamation-triangle"></i> required asset: [`fasta`](available_assets.md#fasta)
<i class="fas fa-tools"></i> required software: [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

```
refgenie build test/bowtie2_index
```

### bismark_bt1_index and bismark_bt2_index

<i class="fas fa-file-import"></i> required files: *none*
<i class="fas fa-sliders-h"></i> required parameters: *none*
<i class="fas fa-exclamation-triangle"></i> required asset: [`fasta`](available_assets.md#fasta)
<i class="fas fa-tools"></i> required software: [bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)

```
refgenie build test/bismark_bt1_index
refgenie build test/bismark_bt2_index
```

### bwa_index

<i class="fas fa-file-import"></i> required files: *none*
<i class="fas fa-sliders-h"></i> required parameters: *none*
<i class="fas fa-exclamation-triangle"></i> required asset: [`fasta`](available_assets.md#fasta)
<i class="fas fa-tools"></i> required software: [bwa](http://bio-bwa.sourceforge.net/)

```
refgenie build test/bwa_index
```

### hisat2_index

<i class="fas fa-file-import"></i> required files: *none*
<i class="fas fa-exclamation-triangle"></i> required asset: [`fasta`](available_assets.md#fasta)
<i class="fas fa-tools"></i> required software: [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml)

```
refgenie build test/hisat2_index
```

### epilog_index

<i class="fas fa-file-import"></i> required files: *none*
<i class="fas fa-sliders-h"></i> required parameters: `--params context=CG` (Default)
<i class="fas fa-exclamation-triangle"></i> required asset: [`fasta`](available_assets.md#fasta)
<i class="fas fa-tools"></i> required software: [epilog](https://github.com/databio/epilog)

```
refgenie build test/epilog_index --params context=CG
```

### kallisto_index

<i class="fas fa-file-import"></i> required files: *none*
<i class="fas fa-sliders-h"></i> required parameters: *none*
<i class="fas fa-exclamation-triangle"></i> required asset: [`fasta`](available_assets.md#fasta)
<i class="fas fa-tools"></i> required software: [kallisto](https://pachterlab.github.io/kallisto/)

```
refgenie build test/kallisto_index
```

### salmon_index

<i class="fas fa-file-import"></i> required files: *none*
<i class="fas fa-sliders-h"></i> required parameters: *none*
<i class="fas fa-exclamation-triangle"></i> required asset: [`fasta`](available_assets.md#fasta)
<i class="fas fa-tools"></i> required software: [salmon](https://salmon.readthedocs.io/en/latest/salmon.html)

```
refgenie build test/salmon_index
```

### star_index

<i class="fas fa-file-import"></i> required files: *none*
<i class="fas fa-sliders-h"></i> required parameters: *none*
<i class="fas fa-exclamation-triangle"></i> required asset: [`fasta`](available_assets.md#fasta)
<i class="fas fa-tools"></i> required software: [star](https://github.com/alexdobin/STAR)

```
refgenie build test/star_index
```

### suffixerator_index

<i class="fas fa-file-import"></i> required files: *none*
<i class="fas fa-sliders-h"></i> required parameters: `--params memlimit=8GB` (Default)
<i class="fas fa-exclamation-triangle"></i> required asset: [`fasta`](available_assets.md#fasta)
<i class="fas fa-tools"></i> required software: [GenomeTools](http://genometools.org/)

```
refgenie build test/suffixerator_index --params memlimit=8GB
```

### tallymer_index

<i class="fas fa-file-import"></i> required files: *none*
<i class="fas fa-sliders-h"></i> required parameters: `--params mersize=30 minocc=2` (Default)
<i class="fas fa-exclamation-triangle"></i> required asset: [`fasta`](available_assets.md#fasta)
<i class="fas fa-tools"></i> required software: [GenomeTools](http://genometools.org/)

```
refgenie build test/tallymer_index --params mersize=30 minocc=2
```

### feat_annotation

<i class="fas fa-file-import"></i> required files: *none*
<i class="fas fa-sliders-h"></i> required parameters: *none*
<i class="fas fa-exclamation-triangle"></i> required asset: [`ensembl_gtf`](build.md#ensembl-gtf), [`ensembl_rb`](build.md#ensembl-rb)
<i class="fas fa-tools"></i> required software: *none*

The `feat_annotation` asset includes the following genomic feature annotations: enhancers, promoters, promoter flanking regions, 5' UTR, 3' UTR, exons, and introns.

```
refgenie build test/feat_annotation
```


[^Amemiya2019]: Amemiya HM, Kundaje A, Boyle AP. The ENCODE Blacklist: Identification of Problematic Regions of the Genome. *Sci Rep* 2019;9, 9354. doi:10.1038/s41598-019-45839-z
