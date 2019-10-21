# This dict provides 'asset packages', which specify recipes (commands) to
# build assets. Each package can produce one or more assets, which are encoded
# as relative paths. The package name is often the same as the asset name but
# it need not be.

# These building recipes should make use of arguments that are auto-populated,
# or user-provided. The auto-populated arguments are:
# - {genome}
# - {asset_outfolder} In addition to these, the recipe should refer in the
#   same way, {var}, to any variables required to be provided, which will be
#   provided via the CLI. These should be listed as 'required_inputs' and
#   will be checked for existence before the commands are executed.

DESC = "description"
ASSET_DESC = "asset_description"
ASSETS = "assets"
PTH = "path"
REQ_IN = "required_inputs"
REQ_ASSETS = "required_assets"
CONT = "container"
CMD_LST = "command_list"


asset_build_packages = {
    "fasta": {
        DESC: "Sequences in the FASTA format, indexed FASTA (produced with samtools index) and chromosome sizes file",
        ASSETS: {
            "fasta": "{genome}.fa",
            "fai": "{genome}.fa.fai",
            "chrom_sizes": "{genome}.chrom.sizes"
        },
        REQ_IN: ["fasta"],
        REQ_ASSETS: [],
        CONT: "databio/refgenie",
        CMD_LST: [
            "cp {fasta} {asset_outfolder}/{genome}.fa.gz",
            "gzip -d {asset_outfolder}/{genome}.fa.gz",
            "samtools faidx {asset_outfolder}/{genome}.fa",
            "cut -f 1,2 {asset_outfolder}/{genome}.fa.fai > {asset_outfolder}/{genome}.chrom.sizes",
        ]
    },
    "dbnsfp": {
        DESC: "A database developed for functional prediction and annotation of all potential non-synonymous single-nucleotide variants (nsSNVs) in the human genome (Gencode release 29/Ensembl 94)",
        ASSETS: {
            "dbnsfp": "{genome}_dbNSFP.txt.gz",
            "tabix": "{genome}_dbNSFP.txt.gz.tbi"
        },
        REQ_IN: ["dbnsfp"],
        REQ_ASSETS: [],
        CONT: "databio/refgenie",
        CMD_LST: [
            "cp {dbnsfp} {asset_outfolder}/{genome}.zip",
            "unzip {asset_outfolder}/{genome}.zip -d {asset_outfolder}",
            "gunzip -v {asset_outfolder}/*variant.chr*.gz",
            "head -n1 {asset_outfolder}/dbNSFP*_variant.chr1 > {asset_outfolder}/{genome}_dbNSFP.txt",
            "cat {asset_outfolder}/dbNSFP*variant.chr* | grep -v '#' >> {asset_outfolder}/{genome}_dbNSFP.txt",
            "rm {asset_outfolder}/dbNSFP*_variant.chr*",
            "bgzip -@ 4 {asset_outfolder}/{genome}_dbNSFP.txt",
            "tabix -s 1 -b 2 -e 2 {asset_outfolder}/{genome}_dbNSFP.txt.gz",
            "rm `find {asset_outfolder} -type f -not -path '{asset_outfolder}/_refgenie_build*' -not -path '{asset_outfolder}/hg38_dbNSFP.txt.*'`"
        ]
    },
    "bowtie2_index": {
        DESC: "Genome index for bowtie, produced with bowtie-build",
        ASSETS: {
            "bowtie2_index": "."
        },
        REQ_IN: [],
        REQ_ASSETS: ["fasta.fasta"],
        CONT: "databio/refgenie",
        CMD_LST: [
            "bowtie2-build {fasta} {asset_outfolder}/{genome}",
            ]
    },
    "bwa_index": {
        DESC: "Genome index for Burrows-Wheeler Alignment Tool, produced with bwa index",
        ASSETS: {
            "bwa_index": "."
        },
        REQ_IN: [],
        REQ_ASSETS: ["fasta.fasta"],
        CONT: "databio/refgenie",
        CMD_LST: [
            "ln -sf {fasta} {asset_outfolder}",
            "bwa index {asset_outfolder}/{genome}.fa",
            ] 
    },    
    "hisat2_index": {
        DESC: "Genome index for HISAT2, produced with hisat2-build",
        ASSETS: {
            "hisat2_index": "."
        },
        REQ_IN: [],
        REQ_ASSETS: ["fasta.fasta"],
        CONT: "databio/refgenie",
        CMD_LST: [
            "hisat2-build {fasta} {asset_outfolder}/{genome}"
            ] 
    },
    "bismark_bt2_index": {
        DESC: "Genome index for Bisulfite-Seq applications, produced by bismark_genome_preparation using bowtie2",
        REQ_IN: [],
        REQ_ASSETS: ["fasta.fasta"],
        CONT: "databio/refgenie",
        ASSETS: {
            "bismark_bt2_index": "."
        },
        CMD_LST: [
            "ln -sf {fasta} {asset_outfolder}",
            "bismark_genome_preparation --bowtie2 {asset_outfolder}"
            ] 
    },
    "bismark_bt1_index": {
        DESC: "Genome index for Bisulfite-Seq applications, produced by bismark_genome_preparation using bowtie1",
        REQ_IN: [],
        REQ_ASSETS: ["fasta.fasta"],
        CONT: "databio/refgenie",
        ASSETS: {
            "bismark_bt1_index": "."
        },
        CMD_LST: [
            "ln -sf {fasta} {asset_outfolder}",
            "bismark_genome_preparation {asset_outfolder}"
            ] 
    },  
    "kallisto_index": {
        DESC: "Genome index for kallisto, produced with kallisto index",
        REQ_IN: [],
        REQ_ASSETS: ["fasta.fasta"],
        CONT: "databio/refgenie",
        ASSETS: {
            "kallisto_index": "."
        },
        CMD_LST: [
            "kallisto index -i {asset_outfolder}/{genome}_kallisto_index.idx {fasta}"
            ] 
    },
    "salmon_index": {
        DESC: "Transcriptome index for salmon, produced with salmon index",
        REQ_IN: [],
        REQ_ASSETS: ["fasta.fasta"],
        CONT: "combinelab/salmon",
        ASSETS: {
            "salmon_index": "."
        },
        CMD_LST: [
            "salmon index -k 31 -i {asset_outfolder} -t {fasta}"
            ] 
    },
    "epilog_index": {
        DESC: "Genome index for CpG sites, produced by the epilog DNA methylation caller",
        REQ_IN: ["context"],
        REQ_ASSETS: ["fasta.fasta"],
        CONT: "databio/refgenie",
        ASSETS: {
            "epilog_index": "."
        },
        CMD_LST: [
            "epilog index -i {fasta} -o {asset_outfolder}/{genome}_{context}.tsv -s {context} -t"
            ] 
    },
    "star_index": {
        DESC: "Genome index for STAR RNA-seq aligner, produced with STAR --runMode genomeGenerate",
        REQ_IN: [],
        REQ_ASSETS: ["fasta.fasta"],
        CONT: "databio/refgenie",
        ASSETS: {
            "star_index": "."
        },
        CMD_LST: [
            "mkdir -p {asset_outfolder}",
            "STAR --runThreadN 16 --runMode genomeGenerate --genomeDir {asset_outfolder} --genomeFastaFiles {fasta}"
            ]
    },
    "gencode_gtf": {
        DESC: "GTF annotation asset which provides access to all annotated transcripts which make up an Ensembl gene set.",
        REQ_IN: ["gencode_gtf"],
        REQ_ASSETS: [],
        CONT: "databio/refgenie",
        ASSETS: {
            "gencode_gtf": "{genome}.gtf.gz"
        },
        CMD_LST: [
            "cp {gencode_gtf} {asset_outfolder}/{genome}.gtf.gz"
            ] 
    },
    "ensembl_gtf": {
        DESC: "Ensembl GTF, TSS, and gene body annotation",
        REQ_IN: ["ensembl_gtf"],
        REQ_ASSETS: [],
        CONT: "databio/refgenie",
        ASSETS: {
            "ensembl_gtf": "{genome}.gtf.gz",
            "ensembl_tss": "{genome}_ensembl_TSS.bed",
            "ensembl_gene_body": "{genome}_ensembl_gene_body.bed",
        },
        CMD_LST: [
            "cp {ensembl_gtf} {asset_outfolder}/{genome}.gtf.gz",
            "gzip -dc {asset_outfolder}/{genome}.gtf.gz | grep 'exon_number \"1\";' | sed 's/^/chr/' | awk -v OFS='\t' '{{print $1, $4, $5, $20, $14, $7}}' | sed 's/\";//g' | sed 's/\"//g' | awk '{{if($6==\"+\"){{print $1\"\t\"$2+20\"\t\"$2+120\"\t\"$4\"\t\"$5\"\t\"$6}}else{{print $1\"\t\"$3-120\"\t\"$3-20\"\t\"$4\"\t\"$5\"\t\"$6}}}}' | LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_ensembl_TSS.bed",
            "gzip -dc {asset_outfolder}/{genome}.gtf.gz | awk '$3 == \"gene\"' | sed 's/^/chr/' | awk -v OFS='\t' '{{print $1,$4,$5,$14,$6,$7}}' | sed 's/\";//g' | sed 's/\"//g' | awk '$4!=\"Metazoa_SRP\"' | awk '$4!=\"U3\"' | awk '$4!=\"7SK\"'  | awk '($3-$2)>200' | awk '{{if($6==\"+\"){{print $1\"\t\"$2+500\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}}else{{print $1\"\t\"$2\"\t\"$3-500\"\t\"$4\"\t\"$5\"\t\"$6}}}}' | awk '$3>$2' | LC_COLLATE=C sort -k4 -u > {asset_outfolder}/{genome}_ensembl_gene_body.bed"
            ] 
    },
    "ensembl_rb": {
        DESC: "A regulatory annotation file",
        REQ_IN: ["gff"],
        REQ_ASSETS: [],
        CONT: "databio/refgenie",
        ASSETS: {
            "ensembl_rb": "{genome}.gff.gz"
        },
        CMD_LST: [
            "cp {gff} {asset_outfolder}/{genome}.gff.gz"
            ] 
    },
    "refgene_anno": {
        DESC: "gene, TSS, exon, intron, and premature mRNA annotation files",
        REQ_IN: ["refgene"],
        REQ_ASSETS: [],
        CONT: "databio/refgenie",
        ASSETS: {
            "refgene_anno": "{genome}_refGene.txt.gz",
            "refgene_tss": "{genome}_TSS.bed",
            "refgene_exon": "{genome}_exons.bed",
            "refgene_intron": "{genome}_introns.bed",
            "refgene_pre_mRNA": "{genome}_pre-mRNA.bed",
        },
        CMD_LST: [
            "cp {refgene} {asset_outfolder}/{genome}_refGene.txt.gz",
            "gzip -dc {asset_outfolder}/{genome}_refGene.txt.gz | awk '{{if($4==\"+\"){{print $3\"\t\"$5\"\t\"$5\"\t\"$13\"\t.\t\"$4}}else{{print $3\"\t\"$6\"\t\"$6\"\t\"$13\"\t.\t\"$4}}}}' | LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_TSS.bed",
            "gzip -dc {asset_outfolder}/{genome}_refGene.txt.gz  | awk -v OFS='\t' '$9>1' | awk -v OFS='\t' '{{ n = split($10, a, \",\"); split($11, b, \",\"); for(i=1; i<n; ++i) print $3, a[i], b[i], $13, i, $4 }}' | awk -v OFS='\t' '$6==\"+\" && $5!=1 {{print $0}} $6==\"-\" {{print $0}}' | awk '$4!=prev4 && prev6==\"-\" {{prev4=$4; prev6=$6; delete line[NR-1]; idx-=1}} {{line[++idx]=$0; prev4=$4; prev6=$6}} END {{for (x=1; x<=idx; x++) print line[x]}}' | LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_exons.bed",
            "gzip -dc {asset_outfolder}/{genome}_refGene.txt.gz  | awk -v OFS='\t' '$9>1' | awk -F'\t' '{{ exonCount=int($9);split($10,exonStarts,\"[,]\"); split($11,exonEnds,\"[,]\"); for(i=1;i<exonCount;i++) {{printf(\"%s\\t%s\\t%s\\t%s\\t%d\\t%s\\n\",$3,exonEnds[i],exonStarts[i+1],$13,($3==\"+\"?i:exonCount-i),$4);}}}}' | LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_introns.bed",
            "gzip -dc {asset_outfolder}/{genome}_refGene.txt.gz  | grep 'cmpl' | awk  '{{print $3\"\t\"$5\"\t\"$6\"\t\"$13\"\t.\t\"$4}}' | LC_COLLATE=C sort -k1,1 -k2,2n -u >  {asset_outfolder}/{genome}_pre-mRNA.bed"
            ]
    },
    "feat_annotation": {
        DESC: "Combined genomic feature annotation created using an Ensembl GTF annotation asset and an Ensembl regulatory build annotation asset",
        ASSETS: {
            "feat_annotation": "{genome}_annotations.bed.gz",
        },
        REQ_IN: [],
        REQ_ASSETS: ["ensembl_gtf.ensembl_gtf", "ensembl_rb.ensembl_rb"],
        CONT: "databio/refgenie",
        CMD_LST: [
            "gzip -dc {ensembl_gtf} | awk '$3==\"exon\"' | grep -v 'pseudogene' | awk -v OFS='\t' '{{print \"chr\"$1, $4-1, $5, \"Exon\", $6, $7}}' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -u > {asset_outfolder}/{genome}_exons.bed",
            "gzip -dc {ensembl_gtf} | awk '$3==\"exon\"' | grep -v 'pseudogene' | awk -v OFS='\t' '{{ split($20, a, \"\\\"\"); print \"chr\"$1, $4-1, $5, a[2], $6, $7}}' | env LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -u | awk 'seen[$4]++ && seen[$4] > 1' | env LC_COLLATE=C sort -k1,1 -k2,2n -k3,3nr | env LC_COLLATE=C sort -k1,1 -k2,2n -u | env LC_COLLATE=C sort -k1,1 -k3,3n -u | awk -v OFS='\t' '{{if($4==prev4){{new2=prev3+1;}} {{prev4=$4; prev3=$3; print $1, new2, $2-1, \"Intron\", $5, $6}}}}' | awk -F'\t' '$2' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_introns.bed",
            "gzip -dc {ensembl_gtf} | awk '$3==\"three_prime_utr\"' | grep -v 'pseudogene' | awk -v OFS='\t' '{{print \"chr\"$1, $4-1, $5, \"3\'\\\'\' UTR\", $6, $7}}' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_3utr.bed",
            "gzip -dc {ensembl_gtf}| awk '$3==\"five_prime_utr\"' | grep -v 'pseudogene' | awk -v OFS='\t' '{{print \"chr\"$1, $4-1, $5, \"5\'\\\'\' UTR\", $6, $7}}' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_5utr.bed",
            "gzip -dc {ensembl_rb} | awk '$3==\"promoter\"' | awk -v OFS='\t' '{{print \"chr\"$1, $4, $5, \"Promoter\", $6, $7}}' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -u > {asset_outfolder}/{genome}_promoter.bed",
            "gzip -dc {ensembl_rb} | awk '$3==\"promoter_flanking_region\"' | awk -v OFS='\t' '{{print \"chr\"$1, $4, $5, \"Promoter Flanking Region\", $6, $7}}' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -u > {asset_outfolder}/{genome}_promoter_flanking.bed",
            "gzip -dc {ensembl_rb} | awk '$3==\"enhancer\"' | awk -v OFS='\t' '{{print \"chr\"$1, $4, $5, \"Enhancer\", $6, $7}}' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -u > {asset_outfolder}/{genome}_enhancer.bed",
            "cat {asset_outfolder}/{genome}_enhancer.bed {asset_outfolder}/{genome}_promoter.bed {asset_outfolder}/{genome}_promoter_flanking.bed {asset_outfolder}/{genome}_5utr.bed {asset_outfolder}/{genome}_3utr.bed {asset_outfolder}/{genome}_exons.bed {asset_outfolder}/{genome}_introns.bed | awk -F'\t' '!seen[$1, $2, $3]++' | env LC_COLLATE=C sort -k4d -k1.4,1V -k2,2n -s > {asset_outfolder}/{genome}_annotations.bed",
            "rm -f {asset_outfolder}/{genome}_enhancer.bed {asset_outfolder}/{genome}_promoter.bed {asset_outfolder}/{genome}_promoter_flanking.bed {asset_outfolder}/{genome}_5utr.bed {asset_outfolder}/{genome}_3utr.bed {asset_outfolder}/{genome}_exons.bed {asset_outfolder}/{genome}_introns.bed",
            "gzip {asset_outfolder}/{genome}_annotations.bed"
            ]
    }
}


