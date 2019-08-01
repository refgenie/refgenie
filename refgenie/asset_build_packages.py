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
        DESC: "Given a gzipped fasta file, produces fasta, fai, and chrom_sizes assets",
        ASSETS: {
            "fasta": {
                PTH: "fasta/{genome}.fa",
                ASSET_DESC: "Sequences in the FASTA format"
            },
            "fai": {
                PTH: "fasta/{genome}.fa.fai",
                ASSET_DESC: "Indexed fasta file, produced with samtools faidx"
            },
            "chrom_sizes": {
                PTH: "fasta/{genome}.chrom.sizes",
                ASSET_DESC: "Chromosome sizes file"
            }
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
    "bowtie2_index": {
        ASSETS: {
            "bowtie2_index": {
                PTH: "bowtie2_index",
                ASSET_DESC: "Genome index for bowtie, produced with bowtie-build"
            }
        },
        REQ_IN: [],
        REQ_ASSETS: ["fasta"],
        CONT: "databio/refgenie",
        CMD_LST: [
            "bowtie2-build {asset_outfolder}/../fasta/{genome}.fa {asset_outfolder}/{genome}",
            ] 
    },
    "bwa_index": {
        ASSETS: {
            "bwa_index": {
                PTH: "bwa_index",
                ASSET_DESC: "Genome index for Burrows-Wheeler Alignment Tool, produced with bwa index"
            }
        },
        REQ_IN: [],
        REQ_ASSETS: ["fasta"],
        CONT: "databio/refgenie",
        CMD_LST: [
            "ln -sf ../fasta/{genome}.fa {asset_outfolder}",
            "bwa index {asset_outfolder}/{genome}.fa",
            ] 
    },    
    "hisat2_index": {
        ASSETS: {
            "hisat2_index": {
                PTH: "hisat2_index",
                ASSET_DESC: "Genome index for HISAT2, produced with hisat2-build"
            }
        },     
        REQ_IN: [],
        REQ_ASSETS: ["fasta"],
        CONT: "databio/refgenie",
        CMD_LST: [
            "hisat2-build {asset_outfolder}/../fasta/{genome}.fa {asset_outfolder}/{genome}"
            ] 
    },
    "bismark_bt2_index": {
        DESC: "The fasta asset must be built first for this to work.",
        REQ_IN: [],
        REQ_ASSETS: ["fasta"],
        CONT: "databio/refgenie",
        ASSETS: {
            "bismark_bt2_index": {
                PTH: "bismark_bt2_index",
                ASSET_DESC: "Genome index for Bisulfite-Seq applications, produced by bismark_genome_preparation using bowtie2"
            }
        },       
        CMD_LST: [
            "ln -sf ../fasta/{genome}.fa {asset_outfolder}",
            "bismark_genome_preparation --bowtie2 {asset_outfolder}"
            ] 
    },
    "bismark_bt1_index": {
        DESC: "The fasta asset must be built first for this to work.",
        REQ_IN: [],
        REQ_ASSETS: ["fasta"],
        CONT: "databio/refgenie",
        ASSETS: {
            "bismark_bt1_index": {
                PTH: "bismark_bt1_index",
                ASSET_DESC: "Genome index for Bisulfite-Seq applications, produced by bismark_genome_preparation using bowtie"
            }
        },       
        CMD_LST: [
            "ln -sf ../fasta/{genome}.fa {asset_outfolder}",
            "bismark_genome_preparation {asset_outfolder}"
            ] 
    },  
    "kallisto_index": {
        REQ_IN: [],
        REQ_ASSETS: ["fasta"],
        CONT: "databio/refgenie",
        ASSETS: {
            "kallisto_index": {
                PTH: "kallisto_index",
                ASSET_DESC: "Genome index for kallisto, produced with kallisto index"
            }
        },
        CMD_LST: [
            "kallisto index -i {asset_outfolder}/{genome}_kallisto_index.idx {asset_outfolder}/../fasta/{genome}.fa"
            ] 
    },
    "salmon_index": {
        REQ_IN: [],
        REQ_ASSETS: ["fasta"],
        CONT: "combinelab/salmon",
        ASSETS: {
            "salmon_index": {
                PTH: "salmon_index",
                ASSET_DESC: "Transcriptome index for salmon, produced with salmon index"
            }
        },
        CMD_LST: [
            "salmon index -k 31 -i {asset_outfolder} -t {asset_outfolder}/../fasta/{genome}.fa"
            ] 
    },
    "epilog_index": {
        REQ_IN: ["context"],
        REQ_ASSETS: ["fasta"],
        CONT: "databio/refgenie",
        ASSETS: {
            "epilog_index": {
                PTH: "epilog_index",
                ASSET_DESC: "Genome index for CpG sites, produced by the epilog DNA methylation caller"
            }
        },
        CMD_LST: [
            "epilog index -i {asset_outfolder}/../fasta/{genome}.fa -o {asset_outfolder}/{genome}_{context}.tsv -s {context} -t"
            ] 
    },
    "star_index": {
        REQ_IN: [],
        REQ_ASSETS: ["fasta"],
        CONT: "databio/refgenie",
        ASSETS: {
            "star_index": {
                PTH: "star_index",
                ASSET_DESC: "Genome index for STAR RNA-seq aligner, produced with STAR --runMode genomeGenerate"
            }
        },
        CMD_LST: [
            "mkdir -p {asset_outfolder}",
            "STAR --runThreadN 16 --runMode genomeGenerate --genomeDir {asset_outfolder} --genomeFastaFiles {asset_outfolder}/../fasta/{genome}.fa "
            ]
    },
    "gencode_gtf": {
        DESC: "Given a GTF file (must be downloaded), create a GTF annotation asset.  GTF provides access to all annotated transcripts which make up an Ensembl gene set.",
        REQ_IN: ["gtf"],
        REQ_ASSETS: [],
        CONT: "databio/refgenie",
        ASSETS: {
            "gencode_gtf": {
                PTH: "gencode_gtf/{genome}.gtf.gz",
                ASSET_DESC: "GTF provides access to all annotated transcripts which make up an Ensembl gene set"
            }
        },
        CMD_LST: [
            "cp {gtf} {asset_outfolder}/{genome}.gtf.gz"
            ] 
    },
    "ensembl_gtf": {
        DESC: "Given a Ensembl GTF file (must be downloaded), create Ensembl GTF, TSS, and gene body annotation assets.",
        REQ_IN: ["gtf"],
        REQ_ASSETS: [],
        CONT: "databio/refgenie",
        ASSETS: {
            # TODO: add asset descriptions
            "ensembl_gtf": {
                PTH: "ensembl_gtf/{genome}.gtf.gz",
                ASSET_DESC: ""
            },
            "ensembl_tss": {
                PTH: "ensembl_gtf/{genome}_ensembl_TSS.bed",
                ASSET_DESC: ""
            },
            "ensembl_gene_body": {
                PTH: "ensembl_gtf/{genome}_ensembl_gene_body.bed",
                ASSET_DESC: ""
            }
        },
        CMD_LST: [
            "cp {gtf} {asset_outfolder}/{genome}.gtf.gz",
            "gzip -dc {asset_outfolder}/../ensembl_gtf/{genome}.gtf.gz | grep 'exon_number \"1\";' | sed 's/^/chr/' | awk -v OFS='\t' '{{print $1, $4, $5, $20, $14, $7}}' | sed 's/\";//g' | sed 's/\"//g' | awk '{{if($6==\"+\"){{print $1\"\t\"$2+20\"\t\"$2+120\"\t\"$4\"\t\"$5\"\t\"$6}}else{{print $1\"\t\"$3-120\"\t\"$3-20\"\t\"$4\"\t\"$5\"\t\"$6}}}}' | LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_ensembl_TSS.bed",
            "gzip -dc {asset_outfolder}/../ensembl_gtf/{genome}.gtf.gz | awk '$3 == \"gene\"' | sed 's/^/chr/' | awk -v OFS='\t' '{{print $1,$4,$5,$14,$6,$7}}' | sed 's/\";//g' | sed 's/\"//g' | awk '$4!=\"Metazoa_SRP\"' | awk '$4!=\"U3\"' | awk '$4!=\"7SK\"'  | awk '($3-$2)>200' | awk '{{if($6==\"+\"){{print $1\"\t\"$2+500\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}}else{{print $1\"\t\"$2\"\t\"$3-500\"\t\"$4\"\t\"$5\"\t\"$6}}}}' | awk '$3>$2' | LC_COLLATE=C sort -k4 -u > {asset_outfolder}/{genome}_ensembl_gene_body.bed"
            ] 
    },
    "ensembl_rb": {
        DESC: "Given a GFF regulatory build file (must be downloaded), create a regulatory annotation asset.",
        REQ_IN: ["gff"],
        REQ_ASSETS: [],
        CONT: "databio/refgenie",
        ASSETS: {
            # TODO: add asset descriptions
            "ensembl_rb": {
                PTH: "ensembl_rb/{genome}.gff.gz",
                ASSET_DESC: ""
            }
        },
        CMD_LST: [
            "cp {gff} {asset_outfolder}/{genome}.gff.gz"
            ] 
    },
    "refgene_anno": {
        DESC: "Given a refGene file (must be downloaded), create gene, TSS, exon, intron, and premature mRNA annotation assets.",
        REQ_IN: ["refgene"],
        REQ_ASSETS: [],
        CONT: "databio/refgenie",
        ASSETS: {
            # TODO: add asset descriptions
            "refgene_anno": {
                PTH: "refgene_anno/{genome}_refGene.txt.gz",
                ASSET_DESC: ""
            },
            "refgene_tss": {
                PTH: "refgene_anno/{genome}_TSS.bed",
                ASSET_DESC: ""
            },
            "refgene_exon": {
                PTH: "refgene_anno/{genome}_exons.bed",
                ASSET_DESC: ""
            },
            "refgene_intron": {
                PTH: "refgene_anno/{genome}_introns.bed",
                ASSET_DESC: ""
            },
            "refgene_pre_mRNA": {
                PTH: "refgene_anno/{genome}_pre-mRNA.bed",
                ASSET_DESC: ""
            }
        },
        CMD_LST: [
            "cp {refgene} {asset_outfolder}/{genome}_refGene.txt.gz",
            "gzip -dc {asset_outfolder}/../refgene_anno/{genome}_refGene.txt.gz | awk '{{if($4==\"+\"){{print $3\"\t\"$5\"\t\"$5\"\t\"$13\"\t.\t\"$4}}else{{print $3\"\t\"$6\"\t\"$6\"\t\"$13\"\t.\t\"$4}}}}' | LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_TSS.bed",
            "gzip -dc {asset_outfolder}/../refgene_anno/{genome}_refGene.txt.gz | awk -v OFS='\t' '$9>1' | awk -v OFS='\t' '{{ n = split($10, a, \",\"); split($11, b, \",\"); for(i=1; i<n; ++i) print $3, a[i], b[i], $13, i, $4 }}' | awk -v OFS='\t' '$6==\"+\" && $5!=1 {{print $0}} $6==\"-\" {{print $0}}' | awk '$4!=prev4 && prev6==\"-\" {{prev4=$4; prev6=$6; delete line[NR-1]; idx-=1}} {{line[++idx]=$0; prev4=$4; prev6=$6}} END {{for (x=1; x<=idx; x++) print line[x]}}' | LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_exons.bed",
            "gzip -dc {asset_outfolder}/../refgene_anno/{genome}_refGene.txt.gz | awk -v OFS='\t' '$9>1' | awk -F'\t' '{{ exonCount=int($9);split($10,exonStarts,\"[,]\"); split($11,exonEnds,\"[,]\"); for(i=1;i<exonCount;i++) {{printf(\"%s\\t%s\\t%s\\t%s\\t%d\\t%s\\n\",$3,exonEnds[i],exonStarts[i+1],$13,($3==\"+\"?i:exonCount-i),$4);}}}}' | LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_introns.bed",
            "gzip -dc {asset_outfolder}/../refgene_anno/{genome}_refGene.txt.gz | grep 'cmpl' | awk  '{{print $3\"\t\"$5\"\t\"$6\"\t\"$13\"\t.\t\"$4}}' | LC_COLLATE=C sort -k1,1 -k2,2n -u >  {asset_outfolder}/{genome}_pre-mRNA.bed"
            ]
    },
    "feat_annotation": {
        DESC: "Using a Ensembl GTF annotation asset and an Ensembl regulatory build annotation asset, create a combined genomic feature annotation asset.",
        ASSETS: {
            # TODO: add asset descriptions
            "feat_annotation": {
                PTH: "feat_annotation/{genome}_annotations.bed.gz",
                ASSET_DESC: ""
            }
        },
        REQ_IN: [],
        REQ_ASSETS: ["ensembl_gtf", "ensembl_rb"],
        CONT: "databio/refgenie",
        CMD_LST: [
            "gzip -dc {asset_outfolder}/../ensembl_gtf/{genome}.gtf.gz | awk '$3==\"exon\"' | grep -v 'pseudogene' | awk -v OFS='\t' '{{print \"chr\"$1, $4-1, $5, \"Exon\", $6, $7}}' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -u > {asset_outfolder}/{genome}_exons.bed",
            "gzip -dc {asset_outfolder}/../ensembl_gtf/{genome}.gtf.gz | awk '$3==\"exon\"' | grep -v 'pseudogene' | awk -v OFS='\t' '{{ split($20, a, \"\\\"\"); print \"chr\"$1, $4-1, $5, a[2], $6, $7}}' | env LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -u | awk 'seen[$4]++ && seen[$4] > 1' | env LC_COLLATE=C sort -k1,1 -k2,2n -k3,3nr | env LC_COLLATE=C sort -k1,1 -k2,2n -u | env LC_COLLATE=C sort -k1,1 -k3,3n -u | awk -v OFS='\t' '{{if($4==prev4){{new2=prev3+1;}} {{prev4=$4; prev3=$3; print $1, new2, $2-1, \"Intron\", $5, $6}}}}' | awk -F'\t' '$2' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_introns.bed",
            "gzip -dc {asset_outfolder}/../ensembl_gtf/{genome}.gtf.gz | awk '$3==\"three_prime_utr\"' | grep -v 'pseudogene' | awk -v OFS='\t' '{{print \"chr\"$1, $4-1, $5, \"3\'\\\'\' UTR\", $6, $7}}' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_3utr.bed",
            "gzip -dc {asset_outfolder}/../ensembl_gtf/{genome}.gtf.gz | awk '$3==\"five_prime_utr\"' | grep -v 'pseudogene' | awk -v OFS='\t' '{{print \"chr\"$1, $4-1, $5, \"5\'\\\'\' UTR\", $6, $7}}' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_5utr.bed",
            "gzip -dc {asset_outfolder}/../ensembl_rb/{genome}.gff.gz | awk '$3==\"promoter\"' | awk -v OFS='\t' '{{print \"chr\"$1, $4, $5, \"Promoter\", $6, $7}}' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -u > {asset_outfolder}/{genome}_promoter.bed",
            "gzip -dc {asset_outfolder}/../ensembl_rb/{genome}.gff.gz | awk '$3==\"promoter_flanking_region\"' | awk -v OFS='\t' '{{print \"chr\"$1, $4, $5, \"Promoter Flanking Region\", $6, $7}}' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -u > {asset_outfolder}/{genome}_promoter_flanking.bed",
            "gzip -dc {asset_outfolder}/../ensembl_rb/{genome}.gff.gz | awk '$3==\"enhancer\"' | awk -v OFS='\t' '{{print \"chr\"$1, $4, $5, \"Enhancer\", $6, $7}}' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -u > {asset_outfolder}/{genome}_enhancer.bed",
            "cat {asset_outfolder}/{genome}_enhancer.bed {asset_outfolder}/{genome}_promoter.bed {asset_outfolder}/{genome}_promoter_flanking.bed {asset_outfolder}/{genome}_5utr.bed {asset_outfolder}/{genome}_3utr.bed {asset_outfolder}/{genome}_exons.bed {asset_outfolder}/{genome}_introns.bed | awk -F'\t' '!seen[$1, $2, $3]++' | env LC_COLLATE=C sort -k4d -k1.4,1V -k2,2n -s > {asset_outfolder}/{genome}_annotations.bed",
            "rm -f {asset_outfolder}/{genome}_enhancer.bed {asset_outfolder}/{genome}_promoter.bed {asset_outfolder}/{genome}_promoter_flanking.bed {asset_outfolder}/{genome}_5utr.bed {asset_outfolder}/{genome}_3utr.bed {asset_outfolder}/{genome}_exons.bed {asset_outfolder}/{genome}_introns.bed",
            "gzip {asset_outfolder}/{genome}_annotations.bed"
            ]
    }
}


