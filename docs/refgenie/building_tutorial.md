# Building and serving tutorial

This tutorial follows an asset from start to finish: set up refgenie, register a genome, build assets, stage them, serve them, and push them to cloud storage. It uses a tiny genome (about 16 KB), so each step runs in seconds.

Refgenie handles an asset in three steps:

```
build  →  stage  →  push
```

| Step | What happens | Result |
|---|---|---|
| **build** | Run a recipe | Asset files in the genome folder, and a database record |
| **stage** | Make a built asset ready to serve | A symlink (file mode) or `.tgz` (archive mode) in the stage folder |
| **push** | Upload staged assets to cloud storage | Files in the cloud, and a record of what was pushed |

You choose what to build, what to stage, and what to push.

## Prerequisites

- Python 3.11 or newer
- refgenie (`pip install "refgenie>=1.0.0a2"`; see [Install](install.md))
- Docker or Apptainer, to run each recipe's tools in its container
- For section 6: `pip install "refgenie[server]>=1.0.0a2"`
- For section 7: `pip install "refgenie[snakemake]>=1.0.0a2"`
- For section 8: a cloud command-line tool, such as `aws`

## 1. Set up refgenie

Create a configuration with its own genome folder and stage folder:

```bash
mkdir -p /tmp/rg_tutorial/{genomes,staged}

refgenie init \
  --genome-folder /tmp/rg_tutorial/genomes \
  --genome-stage-folder /tmp/rg_tutorial/staged
```

A new install has no recipes. Add the refgenie [data channel](data_channels.md) to get them:

```bash
refgenie data-channel add refgenie https https://refgenie.github.io/refgenie-registry/index.yaml
refgenie data-channel sync refgenie --exists-ok
```

Check the result with `refgenie config get`.

## 2. Register a genome

Download rCRSd, a doubled copy of the human mitochondrial reference sequence:

```bash
curl -L https://raw.githubusercontent.com/refgenie/refgenie/master/tests/data/rCRSd.fa \
  -o /tmp/rCRSd.fa
```

Register it:

```bash
refgenie genome init \
  --fasta /tmp/rCRSd.fa \
  --name rCRSd \
  --description "Revised Cambridge Reference Sequence (mitochondrial)" \
  --species "Homo sapiens"
```

Refgenie computes a digest for each sequence and stores the sequences in its RefgetStore, computes a digest for the whole genome, and sets the alias `rCRSd` to point to it. Because a `fasta` recipe is registered, it also builds the genome's `fasta` asset (FASTA file, index, and chromosome sizes). Add `--no-build` to skip that.

```bash
refgenie genome list
refgenie seek rCRSd/fasta              # FASTA file
refgenie seek rCRSd/fasta.fai          # FASTA index
refgenie seek rCRSd/fasta.chrom_sizes  # chromosome sizes
```

See the [Genome initialization tutorial](genome_tutorial.md) for other ways to register a genome.

## 3. Asset classes and recipes

An asset class says what an asset *contains* (its seek keys and serving modes). A recipe says how to *build* one. See [Asset classes and recipes](flexible_asset_types.md).

```bash
refgenie asset-class list
refgenie recipe list
refgenie recipe requirements bowtie2_index
```

The last command shows that `bowtie2_index` needs one input asset, a `fasta`, and an optional `threads` parameter. Recipes take three kinds of input:

| Input | Flag | Example |
|---|---|---|
| asset | `--assets` | `--assets fasta=rCRSd/fasta` |
| file | `--files` | `--files gtf=/path/to/genes.gtf` |
| parameter | `--params` | `--params threads=4` |

## 4. Build an asset

Build a bowtie2 index. Refgenie uses the genome's `fasta` asset as input, and `--container` runs `bowtie2-build` in the recipe's container image:

```bash
refgenie build rCRSd/bowtie2_index --container
```

To name the input asset yourself:

```bash
refgenie build rCRSd/bowtie2_index --container --assets fasta=rCRSd/fasta:default
```

Add `--stage` to build and stage in one step:

```bash
refgenie build rCRSd/bwa_index --container --stage
```

`refgenie list` shows every asset you have.

## 5. Stage assets

Staging prepares a built asset in the stage folder, in the form its asset class's [serving modes](serving_modes.md) ask for:

| Mode | What goes in the stage folder | Use |
|---|---|---|
| `file` | A folder symlink to the built asset | Each file has its own download address |
| `archive` | A `.tgz` archive | One download, for aligner indexes and other bundles |
| `none` | Nothing | Metadata only; users build it themselves |

The `fasta` class uses both `file` and `archive`; `bowtie2_index` uses `archive`.

```bash
refgenie stage add rCRSd/fasta
refgenie stage add rCRSd/bowtie2_index
ls -la /tmp/rg_tutorial/staged/
refgenie stage list
```

To unstage, run `refgenie stage remove rCRSd/fasta`. The built asset stays.

## 6. Serve assets

```bash
refgenie serve              # http://localhost:8000
refgenie serve --port 9000
```

Some of the endpoints:

| Endpoint | What it returns |
|---|---|
| `GET /v4/assets` | All assets |
| `GET /v4/assets/{digest}/files` | The files of a file-mode asset |
| `GET /v4/assets/{digest}/files/{path}` | One file |
| `GET /v4/archives/{digest}/download` | The `.tgz` of an archive-mode asset |
| `GET /ga4gh/drs/v1/objects/{id}` | A GA4GH DRS object |

Get the address of a file without downloading it:

```bash
refgenie seekr rCRSd/fasta -s http://localhost:8000
```

On another computer, subscribe to your server and pull as usual:

```bash
refgenie subscribe -s http://yourserver.example.com
refgenie pull rCRSd/fasta
```

A pull downloads and unpacks the `.tgz` for archive-mode assets, downloads files one by one for file-mode assets, and fails with a message to build locally for `none`-mode assets.

## 7. Build many assets with PEP and Snakemake

To build many assets for many genomes, list them in a [PEP](https://pep.databio.org/) and let refgenie write a [Snakemake](https://snakemake.readthedocs.io/) workflow that builds them in parallel, parents first.

Make a project folder:

```
/tmp/rg_pep/
    config.yaml        # Snakemake config; may be empty
    pep/
        config.yaml
        samples.csv
```

`pep/config.yaml`:

```yaml
pep_version: 2.1.0
sample_table: samples.csv

sample_modifiers:
  derive:
    attributes: [fasta_file_path]
    sources:
      local: /tmp/{genome_name}.fa
```

`pep/samples.csv` has one row per asset to build. `sample_name` is the genome name, and `asset_group_name` must match a recipe name:

```csv
sample_name,genome_name,asset_group_name,fasta_file_path
rCRSd,rCRSd,fasta,local
rCRSd,rCRSd,bowtie2_index,local
rCRSd,rCRSd,bwa_index,local
```

Generate the Snakefile and run it from the project folder. The workflow requires `REFGENIE_INPUTS` to be set:

```bash
cd /tmp/rg_pep
echo "{}" > config.yaml
refgenie generate snakefile --output-path Snakefile
export REFGENIE_INPUTS=/tmp
snakemake --jobs 4
```

The Snakefile has a rule that registers each genome from its FASTA file, and one build rule per recipe that builds and stages the asset. It is written from the recipes registered when you generate it, so generate it again after you add recipes.

On a SLURM cluster, use a [Snakemake profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) with `executor: slurm`:

```bash
snakemake --jobs unlimited --workflow-profile /path/to/snakemake_profile_dir
```

By default the generated build rules run `refgenie build` without `--container`, so each recipe's tools must be on your `PATH`, for example through [bulker](https://bulker.io). To run each recipe in its container image instead, add `--container` to `refgenie generate snakefile`.

## 8. Push to cloud storage

Pushing uploads staged assets to cloud storage (S3 and others) so the public can download them. `refgenie push` uploads the assets you staged with `--push-to`, so add the remote first:

```bash
refgenie remote add \
  --type s3 \
  --prefix my-bucket/refgenie \
  --name production \
  --push-command "aws s3 cp {local_path} s3://{prefix}/{relative_path}"
```

Then build and stage with `--push-to`, and push:

```bash
refgenie build rCRSd/bowtie2_index --container --stage --push-to production
refgenie push --dry-run   # show what would run
refgenie push
refgenie remote status    # pushed and unpushed counts
```

See [Configure remote storage](remotes.md) for the push command placeholders and the `folder_sync` strategy, which uploads the whole stage folder in one command.

## Summary

```bash
# 1. Set up
refgenie init --genome-folder /path/to/genomes --genome-stage-folder /path/to/staged
refgenie data-channel add refgenie https https://refgenie.github.io/refgenie-registry/index.yaml
refgenie data-channel sync refgenie --exists-ok

# 2. Register a genome (also builds its fasta asset)
refgenie genome init --fasta rCRSd.fa --name rCRSd

# 3. Add a remote
refgenie remote add --type s3 --prefix bucket/refgenie --name production \
  --push-command "aws s3 cp {local_path} s3://{prefix}/{relative_path}"

# 4. Build, stage, and mark for push
refgenie build rCRSd/bowtie2_index --container --stage --push-to production

# 5. Serve locally, or push to the cloud
refgenie serve
refgenie push
```

| Command | What it does |
|---|---|
| `refgenie init` | Create the configuration and database |
| `refgenie data-channel add` / `sync` | Get asset classes and recipes |
| `refgenie genome init --fasta FILE --name NAME` | Register a genome from a FASTA file |
| `refgenie recipe requirements NAME` | Show a recipe's inputs |
| `refgenie build GENOME/ASSET` | Build an asset |
| `refgenie build GENOME/ASSET --stage` | Build and stage |
| `refgenie list` | List local assets |
| `refgenie seek GENOME/ASSET` | Get a local file path |
| `refgenie stage add` / `remove` / `list` | Manage staged assets |
| `refgenie serve` | Start a server |
| `refgenie seekr GENOME/ASSET -s URL` | Get a file's address on a server |
| `refgenie generate snakefile -o FILE` | Write a Snakemake build workflow |
| `refgenie remote add` / `list` / `status` | Manage remotes |
| `refgenie push` | Upload assets marked for push |
