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


TSS_AWK = "'{if($4==\"+\"){print $3\"\t\"$5\"\t\"$5\"\t\"$4\"\t\"$13}else{print $3\"\t\"$6\"\t\"$6\"\t\"$4\"\t\"$13}}'"
GENE_ANNO_FILE_NAME = "refGene.txt"
GENE_ANNO_NAME = "gene_anno"

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
        "required_inputs": ["gtf"],
        "required_assets": [],
        "container": "databio/refgenie",
        "assets": {
            "gtf_anno": "gtf_anno"
            },
        "command_list": [
            "cp {gtf} {asset_outfolder}/{genome}.gtf",
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
    GENE_ANNO_NAME: {
        "required_inputs": ["annogene"],
        "required_assets": [],
        #"container": "databio/refgenie",
        "assets": {
            GENE_ANNO_NAME: GENE_ANNO_NAME
        },
        "command_list": [
            "cp {{annogene}} {{asset_outfolder}}/{{genome}}_{annsfile}".format(annsfile=GENE_ANNO_FILE_NAME)
        ]
    },
    "TSS_annotation": {
        "required_inputs": [],
        "required_assets": [GENE_ANNO_NAME],
        #"container": "databio/refgenie",
        "command_list": [
            "gunzip -c {{asset_outfolder}}/{{genome}}_{annsfile} | awk {awk_cmd} | LC_COLLATE=C sort -k1,1 -k2,2n -u > {{asset_outfolder}}{{genome}}_TSS.tsv".
                format(annsfile=GENE_ANNO_FILE_NAME, awk_cmd=TSS_AWK)
        ]
    }
}


