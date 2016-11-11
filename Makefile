ROOT=${RESOURCES}genomes/
BUILDER=${CODEBASE}pipelines/pipelines/build_reference.py


hg19-1: $(ROOT)hg19/ $(ROOT)hg19/hg19.2bit $(ROOT)/hg19/hg19.gtf.gz
	twoBitToFa $(ROOT)hg19/hg19.2bit $(ROOT)hg19/hg19.fa
	$${CODEBASE}pipelines/pipelines/build_reference.py -f $(ROOT)/hg19/hg19.fa -a $(ROOT)/hg19/hg19.gtf.gz

$(ROOT)/hg19/:
	mkdir -p $(ROOT)hg19/

$(ROOT)/hg19/hg19.2bit:
	wget -O $(ROOT)/hg19/hg19.2bit http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit

$(ROOT)/hg19/hg19.gtf.gz:
	wget -O $(ROOT)/hg19/hg19.gtf.gz ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz


hg19:
	INPUT=http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
	GTF=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
	${BUILDER} -f ${INPUT} -a ${GTF} -n hg19

# 	# name used internally:
# GENOME_NAME=hg38
# # use the NCBI's official version for sequence alignments without _alt sequences:
# FASTA_SOURCE=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
# # get genes from GENCODE? (optional)
# # (http://www.gencodegenes.or; use the version for the "primary assembly", 
# # which contains chromosomes + scaffolds, but not _alt assemblies -- 
# # this should match the version we used for the genome):
# GTF_SOURCE=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_23/gencode.v23.primary_assembly.annotation.gtf.gz


# ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz
