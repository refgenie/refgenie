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

asset_build_packages = {
    "fasta": {
        "description": "Given a gzipped fasta file, produces fasta, fai, and chrom_sizes assets",
        "assets": {
            "fasta": "fasta/{genome}.fa",
            "fai": "fasta/{genome}.fa.fai",
            "chrom_sizes": "fasta/{genome}.chrom.sizes",
        },
        "required_inputs": ["fasta"],
        "required_assets": [],
        "container": "databio/refgenie",
        "command_list": [
            "cp {fasta} {asset_outfolder}/{genome}.fa.gz",
            "gzip -d {asset_outfolder}/{genome}.fa.gz",
            "samtools faidx {asset_outfolder}/{genome}.fa",
            "cut -f 1,2 {asset_outfolder}/{genome}.fa.fai > {asset_outfolder}/{genome}.chrom.sizes",
        ]
    },
    "bowtie2_index": {
        "assets": {
            "bowtie2_index": "bowtie2_index",
        },
        "required_inputs": [],
        "required_assets": ["fasta"],
        "container": "databio/refgenie",
        "command_list": [
            "bowtie2-build {asset_outfolder}/../fasta/{genome}.fa {asset_outfolder}/{genome}",
            ] 
    },
    "bwa_index": {
        "assets": {
            "bwa_index": "bwa_index",
        },
        "required_inputs": [],
        "required_assets": ["fasta"],
        "container": "databio/refgenie",
        "command_list": [
            "ln -sf ../fasta/{genome}.fa {asset_outfolder}",
            "bwa index {asset_outfolder}/{genome}.fa",
            ] 
    },    
    "hisat2_index": {
        "assets": {
            "hisat2_index": "hisat2_index",
        },     
        "required_inputs": [],
        "required_assets": ["fasta"],
        "container": "databio/refgenie",
        "command_list": [
            "hisat2-build {asset_outfolder}/../fasta/{genome}.fa {asset_outfolder}/{genome}"
            ] 
    },
    "bismark_bt2_index": {
        "description": "The fasta asset must be built first for this to work.",
        "required_inputs": [],
        "required_assets": ["fasta"],
        "container": "databio/refgenie",
        "assets": {
            "bismark_bt2_index": "bismark_bt2_index",
        },       
        "command_list": [
            "ln -sf ../fasta/{genome}.fa {asset_outfolder}",
            "bismark_genome_preparation --bowtie2 {asset_outfolder}"
            ] 
    },
    "bismark_bt1_index": {
        "description": "The fasta asset must be built first for this to work.",
        "required_inputs": [],
        "required_assets": ["fasta"],
        "container": "databio/refgenie",
        "assets": {
            "bismark_bt1_index": "bismark_bt1_index",
        },       
        "command_list": [
            "ln -sf ../fasta/{genome}.fa {asset_outfolder}",
            "bismark_genome_preparation {asset_outfolder}"
            ] 
    },  
    "kallisto_index": {
        "required_inputs": [],
        "required_assets": ["fasta"],
        "container": "databio/refgenie",
        "assets": {
            # "kallisto_index": "{asset_outfolder}/{genome}_kallisto_index.idx"
            "kallisto_index": "kallisto_index"
            },
        "command_list": [
            "kallisto index -i {asset_outfolder}/{genome}_kallisto_index.idx {asset_outfolder}/../fasta/{genome}.fa"
            ] 
    },
    "salmon_index": {
        "required_inputs": [],
        "required_assets": ["fasta"],
        "container": "combinelab/salmon",
        "assets": {
            "salmon_index": "salmon_index"
            },
        "command_list": [
            "salmon index -k 31 -i {asset_outfolder} -t {asset_outfolder}/../fasta/{genome}.fa"
            ] 
    },
    "gtf_anno": {
        "description": "Given a GTF file (must be downloaded), create a GTF annotation asset.  GTF provides access to all annotated transcripts which make up an Ensembl gene set.",
        "required_inputs": ["gtf"],
        "required_assets": [],
        "container": "databio/refgenie",
        "assets": {
            "gtf_anno": "gtf_anno/{genome}.gtf.gz"
            },
        "command_list": [
            "cp {gtf} {asset_outfolder}/{genome}.gtf.gz",
            ] 
    },
    "reg_anno": {
        "description": "Given a GTF regulatory file (must be downloaded), create a regulatory annotation asset.",
        "required_inputs": ["gff"],
        "required_assets": [],
        "container": "databio/refgenie",
        "assets": {
            "reg_anno": "reg_anno/{genome}.gff.gz"
            },
        "command_list": [
            "cp {gff} {asset_outfolder}/{genome}.gff.gz",
            ] 
    },
    "gene_anno": {
        "description": "Given a refGene file (must be downloaded), create a gene annotation asset.",
        "required_inputs": ["refgene"],
        "required_assets": [],
        #"container": "databio/refgenie",
        "assets": {
            "gene_anno": "gene_anno/{genome}_refGene.txt.gz"
            },
        "command_list": [
            "cp {refgene} {asset_outfolder}/{genome}_refGene.txt.gz"
            ]
    },
    "epilog_index": {
        "required_inputs": ["context"],
        "required_assets": ["fasta"],
        "container": "databio/refgenie",
        "assets": {
            "epilog_index": "epilog_index"
            },
        "command_list": [
            "epilog index -i {asset_outfolder}/../fasta/{genome}.fa -o {asset_outfolder}/{genome}_{context}.tsv -s {context} -t"
            ] 
    },
    "star_index": {
        "required_inputs": [],
        "required_assets": ["fasta"],
        "container": "databio/refgenie",
        "assets": {
            "star_index": "star_index"
            },
        "command_list": [
            "mkdir -p {asset_outfolder}",
            "STAR --runThreadN 16 --runMode genomeGenerate --genomeDir {asset_outfolder} --genomeFastaFiles {asset_outfolder}/../fasta/{genome}.fa "
            ]
    },
    "tss_annotation": {
        "description": "Using a gene annotation file, create a TSS annotation asset.",
        "required_inputs": [],
        "required_assets": ["gene_anno"],
        "assets": {
            "tss_annotation": "tss_annotation/{genome}_TSS.bed"
            },
        #"container": "databio/refgenie",
        "command_list": [
            "gzip -dc {asset_outfolder}/../gene_anno/{genome}_refGene.txt.gz | awk '{{if($4==\"+\"){{print $3\"\t\"$5\"\t\"$5\"\t\"$13\"\t.\t\"$4}}else{{print $3\"\t\"$6\"\t\"$6\"\t\"$13\"\t.\t\"$4}}}}' | LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_TSS.bed"
            ]
    },
    "exon_annotation": {
        "description": "Using a gene annotation file, create an exon (no first exons) annotation asset.",
        "required_inputs": [],
        "required_assets": ["gene_anno"],
        "assets": {
            "exon_annotation": "exon_annotation/{genome}_exons.bed"
            },
        #"container": "databio/refgenie",
        "command_list": [
            "gzip -dc {asset_outfolder}/../gene_anno/{genome}_refGene.txt.gz | awk -v OFS='\t' '$9>1' | awk -v OFS='\t' '{{ n = split($10, a, \",\"); split($11, b, \",\"); for(i=1; i<n; ++i) print $3, a[i], b[i], $13, i, $4 }}' | awk -v OFS='\t' '$6==\"+\" && $5!=1 {{print $0}} $6==\"-\" {{print $0}}' | awk '$4!=prev4 && prev6==\"-\" {{prev4=$4; prev6=$6; delete line[NR-1]; idx-=1}} {{line[++idx]=$0; prev4=$4; prev6=$6}} END {{for (x=1; x<=idx; x++) print line[x]}}' | LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_exons.bed"
            ]
    },
    "intron_annotation": {
        "description": "Using a gene annotation file, create an intron annotation asset.",
        "required_inputs": [],
        "required_assets": ["gene_anno"],
        "assets": {
            "intron_annotation": "intron_annotation/{genome}_introns.bed"
            },
        #"container": "databio/refgenie",
        "command_list": [
            "gzip -dc {asset_outfolder}/../gene_anno/{genome}_refGene.txt.gz | awk -v OFS='\t' '$9>1' | awk -F'\t' '{{ exonCount=int($9);split($10,exonStarts,\"[,]\"); split($11,exonEnds,\"[,]\"); for(i=1;i<exonCount;i++) {{printf(\"%s\\t%s\\t%s\\t%s\\t%d\\t%s\\n\",$3,exonEnds[i],exonStarts[i+1],$13,($3==\"+\"?i:exonCount-i),$4);}}}}' | LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_introns.bed"
            ]
    },
    "pre_mRNA_annotation": {
        "description": "Using a gene annotation file, create a premature mRNA annotation asset.",
        "required_inputs": [],
        "required_assets": ["gene_anno"],
        "assets": {
            "pre_mRNA_annotation": "pre_mRNA_annotation/{genome}_pre-mRNA.bed"
            },
        #"container": "databio/refgenie",
        "command_list": [
            "gzip -dc {asset_outfolder}/../gene_anno/{genome}_refGene.txt.gz | grep 'cmpl' | awk  '{{print $3\"\t\"$5\"\t\"$6\"\t\"$13\"\t.\t\"$4}}' | LC_COLLATE=C sort -k1,1 -k2,2n -u >  {asset_outfolder}/{genome}_pre-mRNA.bed"
            ]
    },
    "pi_tss": {
        "description": "Using a GTF annotation file, create a pause index *all* possible TSS annotation asset.",
        "required_inputs": [],
        "required_assets": ["gtf_anno"],
        "assets": {
            "pi_tss": "pi_tss/{genome}_PI_TSS.bed"
            },
        #"container": "databio/refgenie",
        "command_list": [
            "gzip -dc {asset_outfolder}/../gtf_anno/{genome}.gtf.gz | grep 'exon_number \"1\"' | sed 's/^/chr/' | awk -v OFS='\t' '{{print $1, $4, $5, $20, $14, $7}}' | sed 's/\";//g' | sed 's/\"//g' | awk '{{if($6==\"+\"){{print $1\"\t\"$2+20\"\t\"$3+120\"\t\"$4\"\t\"$5\"\t\"$6}}else{{print $1\"\t\"$3-120\"\t\"$3-20\"\t\"$4\"\t\"$5\"\t\"$6}}}}' | LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_PI_TSS.bed"
            ]
    },
    "pi_body": {
        "description": "Using a GTF annotation file, create a pause index *all* possible gene body annotation asset.",
        "required_inputs": [],
        "required_assets": ["gtf_anno"],
        "assets": {
            "pi_body": "pi_tss/{genome}_PI_gene_body.bed"
            },
        #"container": "databio/refgenie",
        "command_list": [
            "gzip -dc {asset_outfolder}/../gtf_anno/{genome}.gtf.gz | awk '$3 == \"gene\"' | sed 's/^/chr/' | awk -v OFS='\t' '{{print $1,$4,$5,$14,$6,$7}}' | sed 's/\";//g' | sed 's/\"//g' | awk '$4!=\"Metazoa_SRP\"' | awk '$4!=\"U3\"' | awk '$4!=\"7SK\"'  | awk '($3-$2)>200' | awk '{{if($6==\"+\"){{print $1\"\t\"$2+500\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6}}else{{print $1\"\t\"$2\"\t\"$3-500\"\t\"$4\"\t\"$5\"\t\"$6}}}}' | awk '$3>$2' | LC_COLLATE=C sort -k4 -u > {asset_outfolder}/{genome}_PI_gene_body.bed"
            ]
    },
    "feat_annotation": {
        "description": "Using a GTF annotation file and a regulatory annotation asset, create a combined genomic feature annotation asset.",
        "assets": {
            "exons": "feat_annotation/{genome}_exons.bed",
            "introns": "feat_annotation/{genome}_introns.bed",
            "utr5": "feat_annotation/{genome}_5utr.bed",
            "utr3": "feat_annotation/{genome}_3utr.bed",
            "promoter": "feat_annotation/{genome}_promoter.bed",
            "promoter_flanking": "feat_annotation/{genome}_promoter_flanking.bed",
            "enhancer": "feat_annotation/{genome}_enhancer.bed",
            "feat_annotation": "{asset_outfolder}/{genome}_annotations.bed.gz"
            },
        "required_inputs": [],
        "required_assets": ["gtf_anno", "reg_anno"],
        #"container": "databio/refgenie",
        "command_list": [
            "gzip -dc {asset_outfolder}/../gtf_anno/{genome}.gtf.gz | awk '$3==\"exon\"' | grep -v 'pseudogene' | awk -v OFS='\t' '{{print \"chr\"$1, $4-1, $5, \"Exon\", $6, $7}}' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -u > {asset_outfolder}/{genome}_exons.bed",
            "gzip -dc {asset_outfolder}/../gtf_anno/{genome}.gtf.gz | awk '$3==\"exon\"' | grep -v 'pseudogene' | awk -v OFS='\t' '{{ split($20, a, \"\\\"\"); print \"chr\"$1, $4-1, $5, a[2], $6, $7}}' | env LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -u | awk 'seen[$4]++ && seen[$4] > 1' | env LC_COLLATE=C sort -k1,1 -k2,2n -k3,3nr | env LC_COLLATE=C sort -k1,1 -k2,2n -u | env LC_COLLATE=C sort -k1,1 -k3,3n -u | awk -v OFS='\t' '{{if($4==prev4){{new2=prev3+1;}} {{prev4=$4; prev3=$3; print $1, new2, $2-1, \"Intron\", $5, $6}}}}' | awk -F'\t' '$2' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_introns.bed",
            "gzip -dc {asset_outfolder}/../gtf_anno/{genome}.gtf.gz | awk '$3==\"three_prime_utr\"' | grep -v 'pseudogene' | awk -v OFS='\t' '{{print \"chr\"$1, $4-1, $5, \"3\'\\\'\' UTR\", $6, $7}}' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_3utr.bed",
            "gzip -dc {asset_outfolder}/../gtf_anno/{genome}.gtf.gz | awk '$3==\"five_prime_utr\"' | grep -v 'pseudogene' | awk -v OFS='\t' '{{print \"chr\"$1, $4-1, $5, \"5\'\\\'\' UTR\", $6, $7}}' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -u > {asset_outfolder}/{genome}_5utr.bed",
            "gzip -dc {asset_outfolder}/../reg_anno/{genome}.gff.gz | awk '$3==\"promoter\"' | awk -v OFS='\t' '{{print \"chr\"$1, $4, $5, \"Promoter\", $6, $7}}' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -u > {asset_outfolder}/{genome}_promoter.bed",
            "gzip -dc {asset_outfolder}/../reg_anno/{genome}.gff.gz | awk '$3==\"promoter_flanking_region\"' | awk -v OFS='\t' '{{print \"chr\"$1, $4, $5, \"Promoter Flanking Region\", $6, $7}}' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -u > {asset_outfolder}/{genome}_promoter_flanking.bed",
            "gzip -dc {asset_outfolder}/../reg_anno/{genome}.gff.gz | awk '$3==\"enhancer\"' | awk -v OFS='\t' '{{print \"chr\"$1, $4, $5, \"Enhancer\", $6, $7}}' | awk '$2<$3' | env LC_COLLATE=C sort -k1,1 -k2,2n -k3,3n -u > {asset_outfolder}/{genome}_enhancer.bed",
            "cat {asset_outfolder}/{genome}_enhancer.bed {asset_outfolder}/{genome}_promoter.bed {asset_outfolder}/{genome}_promoter_flanking.bed {asset_outfolder}/{genome}_5utr.bed {asset_outfolder}/{genome}_3utr.bed {asset_outfolder}/{genome}_exons.bed {asset_outfolder}/{genome}_introns.bed | awk -F'\t' '!seen[$1, $2, $3]++' | env LC_COLLATE=C sort -k4d -k1.4,1V -k2,2n -s > {asset_outfolder}/{genome}_annotations.bed",
            "gzip {asset_outfolder}/{genome}_annotations.bed"
            ]
    }
}


