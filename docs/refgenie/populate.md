# Make configuration files portable with refgenie populate

`refgenie populate` replaces registry paths such as `refgenie://hg38/fasta` in text with local file paths. `refgenie populater` replaces them with the asset's address on a server instead (an `http://` or `s3://` URL), for cloud or short-lived compute.

This lets a workflow that knows nothing about refgenie still use refgenie-managed files. For example, [Common Workflow Language](https://www.commonwl.org/) (CWL) workflows need every input path before they start, so you fill in the input file with refgenie first.

## Strings

Pipe text through `refgenie populate`:

```console
echo 'bowtie2 -x refgenie://hg38/bowtie2_index -U r1.fq -S eg1.sam' | refgenie populate | sh
```

## Files

Given `config_template.yaml`:

```yaml
config:
  param1: value1
  fasta: "refgenie://hg38/fasta"
  bowtie2_index: "refgenie://hg38/bowtie2_index"
```

run:

```console
refgenie populate --file config_template.yaml > config.yaml
```

to get:

```yaml
config:
  param1: value1
  fasta: "/home/johndoe/genomes/alias/hg38/fasta/default/hg38.fa"
  bowtie2_index: "/home/johndoe/genomes/alias/hg38/bowtie2_index/2.5.4/hg38"
```

`refgenie populater` takes the same input, and `-s` to choose servers.

## Looper

[Looper](https://looper.databio.org/) can run refgenie before it submits each job. Add the hook to your pipeline interface:

```yaml
pre_submit:
  python_functions:
  - refgenie.integrations.looper.populate
```

Looper then replaces `refgenie://` paths anywhere in the pipeline block, and templates can read paths such as `{ refgenie[sample.genome].fasta.chrom_sizes }`. To use a database other than the default, set `var_templates.refgenie_db_config` in the pipeline interface, or `REFGENIE_DB_CONFIG_PATH`.
