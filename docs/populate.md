# Replace refgenie registry paths with asset file paths 

`refgenie populate` and `refgenie populater` commands can be used to replace refgenie registry paths with asset file paths -- **go from `refgenie://hg38/fasta` to `/home/johndoe/genomes/hg38/fasta/default/hg38.fa`**.

# Motivation

In some cases it is desirable to run a refgenie-unaware workflow and benefit from refgenie framework. Is such cases we need to populate input for a workflow run as a pre-processing step, outside of the workflow. For instance, this is the way [Common Workflow Language](https://www.commonwl.org/) (CWL) works; CWL workflows in best practices require knowledge of all input files before the workflow run begins. So, it doesn't make sense to pass a registry path, which is then resolved by refgenie inside the workflow.

# Usage examples

Both `populate` commands can populate refgenie registry paths both in file and string.

## Text intput

To populate an argument to `bowtie2` command with a local path to `bowtie2_index` asset managed by refgenie and run the aligner call:

```console
echo 'bowtie2 -x refgenie://hg38/bowtie2_index -U reads_1.fq -S eg1.sam' | refgenie populate | sh
```

## File input

Example input in `test/config_template.yaml`:
```yaml
config:
  param1: value1
  fasta: refgenie://hg38/fasta
  bowtie2_index: refgenie://hg38/bowtie2_index
```

To populate a bowtie2 index and FASTA file paths in a YAML configuration file of an arbitrary pipelne call:

```console
refgenie populate --file test/config_template.yaml > test/config.yaml
```

Example output in `test/config.yaml`:
```yaml
config:
  param1: value1
  fasta: /home/johndoe/genomes/hg38/fasta/default/hg38.fa
  bowtie2_index: /home/johndoe/genomes/hg38/bowtie2_index/default/hg38
```
