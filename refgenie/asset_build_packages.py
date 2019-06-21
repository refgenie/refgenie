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
    "fastagz": {
        "assets": {
            "fasta": "fasta/{genome}.fa.gz",
            "fai": "fasta/{genome}.fa.gz.fai",
            "chrom_sizes": "fasta/{genome}.chrom.sizes",
        },
        "required_inputs": ["fasta"],
        "command_list": [
            "cp {fasta} {asset_outfolder}/{genome}.fa.gz",
            "samtools faidx {asset_outfolder}/{genome}.fa.gz",
            "cut -f 1,2 {asset_outfolder}/{genome}.fa.gz.fai > {asset_outfolder}/{genome}.chrom.sizes"
        ]
    },
    "bowtie2_index": {
        "assets": {
            "bowtie2_index": "bowtie2_index",
        },       
        "required_inputs": ["fasta"],
        "command_list": [
            "bowtie2-build {fasta} {asset_outfolder}/{genome}"
            ] 
    },
    "hisat2_index": {
        "assets": {
            "hisat2_index": "hisat2_index",
        },       
        "required_inputs": ["fasta"],
        "command_list": [
            "hisat2-build {fasta} {asset_outfolder}/{genome}"
            ] 
    },
    "bismark_bt2_index": {
        "assets": {
            "bismark_bt2_index": "bismark_bt2_index",
        },       
        "required_inputs": ["fasta"],
        "command_list": [
            "ln -sf ../{genome}.fa.gz {asset_outfolder}"
            "bismark_genome_preparation --bowtie2 {fasta} {asset_outfolder}/{genome}"
            ] 
    },
    "bismark_bt1_index": {
        "assets": {
            "bismark_bt1_index": "bismark_bt1_index",
        },       
        "required_inputs": ["fasta"],
        "command_list": [
            "ln -sf ../{genome}.fa.gz {asset_outfolder}"
            "bismark_genome_preparation {fasta} {asset_outfolder}/{genome}"
            ] 
    },  
    "kallisto_index": {
        "required_inputs": ["fasta"],
        "assets": {
            "kallisto_index": "kallisto_index"
            },
        "command_list": [
            "kallisto index -i {fasta} {asset_outfolder}/{genome}_kallisto_index.idx"
            ] 
    },
    "gtf_anno": {
        "required_inputs": ["gtf"],
        "assets": {
            "gtf_anno": "gtf_anno"
            },
        "command_list": [
            "cp {gtf} {asset_outfolder}/{genome}.gtf",
            ] 
    },
    "epilog_index": {
        "required_inputs": ["fasta", "context"],
        "assets": {
            "epilog_index": "epilog_index"
            },
        "command_list": [
            "epilog index -i {fasta} -o {asset_outfolder}/{genome}_{context}.tsv -s {context} -t"
            ] 
    }
}
