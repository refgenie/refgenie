# Make configuration files portable with refgenie populate commands

Use `refgenie populate` to replace registry paths (*e.g.* `refgenie://hg38/fasta`) in text files with asset file paths (*e.g.* `/home/johndoe/genomes/hg38/fasta/default/hg38.fa`). The remote version, `refgenie populatr`, will replace your registry path with a URI, like `s3://path/to/asset.xyz` or `http://path/to/asset.xyz` for use in an ephemeral compute environment. This powerful features allows you to write configuration files and scripts with maximum portability for anything you might need to configure with reference genome paths.

# Motivation

In some cases it is desirable to run a refgenie-unaware workflow and benefit from the refgenie framework. In such cases we need to populate an input for a workflow run as a pre-processing step, outside of the workflow. For instance, this is the way [Common Workflow Language](https://www.commonwl.org/) (CWL) works; CWL workflows in best practices require knowledge of all input files before the workflow run begins. So, it doesn't make sense to pass a registry path, which is then resolved by refgenie inside the workflow.

# Usage examples

Both `populate` commands can populate refgenie registry paths either in a **file** or a **string**.

## String intput

Use a pipe (`|`) to populate an in-line command argument with a local path managed by refgenie:

```console
echo 'bowtie2 -x refgenie://hg38/bowtie2_index -U reads_1.fq -S eg1.sam' | refgenie populate | sh
```

## File input

Example input in `test/config_template.yaml`:
```yaml
config:
  param1: value1
  fasta: "refgenie://hg38/fasta"
  bowtie2_index: "refgenie://hg38/bowtie2_index"
```

To populate a bowtie2 index and FASTA file paths in a YAML configuration file of an arbitrary pipeline call:

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
