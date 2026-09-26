# Building and Serving Tutorial

This tutorial walks through the complete refgenie asset lifecycle: initialize a genome, build assets, stage them for serving, serve them locally, and push them to cloud storage. Each stage builds on the previous one, and you can follow along using a small test genome (~16 KB) so the steps run quickly.

!!! warning "Transitionary period"
    During the current transition period, the package is installed as `refgenie1` and the CLI command is `refgenie1`. When the transition is complete, it will become `refgenie`. All examples below use `refgenie1`.

## Overview: the three-stage lifecycle

Refgenie separates asset management into three distinct stages:

```
build  →  stage  →  push
```

| Stage | What happens | Output |
|-------|-------------|--------|
| **build** | Run a recipe; produce asset files in `genome_folder` | Files on disk, DB record |
| **stage** | Make a built asset servable in `genome_stage_folder` | Symlink (file mode) or `.tgz` (archive mode), DB StagedAsset record |
| **push** | Upload staged assets to cloud storage | RemoteAssetLink records, cloud files |

This separation gives you control: you can build locally, stage selectively, and push only what you want to host publicly.

## Prerequisites

- Python 3.10+
- refgenie1 installed (`pip install refgenie1-*.whl`)
- Snakemake installed (`pip install snakemake`) — only needed for Section 7
- Cloud CLI (e.g., `aws`, `az`, `gsutil`) — only needed for Section 8

## Section 1: Initialize refgenie

Set up refgenie with explicit paths for where to store assets and staged (servable) assets. Using a temporary directory makes this tutorial self-contained:

```bash
mkdir -p /tmp/rg_tutorial/{genomes,staged}

refgenie1 init \
  --genome-folder /tmp/rg_tutorial/genomes \
  --genome-stage-folder /tmp/rg_tutorial/staged
```

Verify the configuration:

```bash
refgenie1 config get
```

Expected output:

```
                                     Refgenie configuration
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━┓
┃ genome_folder                      ┃ version ┃ genome_stage_folder                ┃ servers ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━┩
│ /tmp/rg_tutorial/genomes           │ 1       │ /tmp/rg_tutorial/staged            │ []      │
└────────────────────────────────────┴─────────┴────────────────────────────────────┴─────────┘
```

## Section 2: Initialize a genome

Download the rCRSd FASTA — the revised Cambridge Reference Sequence for human mitochondrial DNA, ~16 KB:

```bash
curl -L https://raw.githubusercontent.com/refgenie/refgenie-build-fasta/main/test_data/rCRSd.fa \
  -o /tmp/rCRSd.fa
```

Register the genome with refgenie. This reads the FASTA, computes sequence digests, builds a RefgetStore, and records the genome in the database:

```bash
refgenie1 genome init \
  --fasta /tmp/rCRSd.fa \
  --name rCRSd \
  --description "Revised Cambridge Reference Sequence (mitochondrial)" \
  --species "Homo sapiens"
```

What happened:

1. Each sequence in the FASTA was digested (SHA-512/24) and stored in the RefgetStore
2. A sequence collection digest was computed from all sequences
3. The genome was registered in the database with that digest
4. The alias `rCRSd` was set, pointing to the genome digest

Verify:

```bash
refgenie1 genome list
```

!!! tip "See genome_tutorial.md for more"
    For deeper coverage of genome initialization including remote initialization from a seqcolapi server, see the [Genome initialization tutorial](genome_tutorial.md).

## Section 3: Asset classes and recipes

Asset classes define what an asset *contains* (seek keys and serving modes). Recipes define *how to build* one. They are separate so the same recipe can produce assets of different classes, and the same class can be built by different recipes.

List asset classes and recipes available locally:

```bash
refgenie1 asset_class list
refgenie1 recipe list
```

View the requirements for a recipe with the `-q` flag:

```bash
refgenie1 recipe show fasta -q
```

This shows the recipe's input asset classes, input files, input parameters, and associated Docker image.

For a recipe with asset dependencies:

```bash
refgenie1 recipe show bowtie2_index -q
```

Expected output (abbreviated):

```
                                     Recipes
┏━━━━━━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┓
┃ Name          ┃ Version ┃ Output asset class ┃ Input asset classes       ┃ Input files ┃
┡━━━━━━━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━┩
│ bowtie2_index │ 0.0.1   │ bowtie2_index      │ • fasta (fasta asset)    │ None        │
└───────────────┴─────────┴────────────────────┴──────────────────────────┴─────────────┘
```

The three input types:

| Input type | CLI flag | Description |
|-----------|---------|-------------|
| `required_assets` | `--assets` | A pre-built refgenie asset (e.g., fasta) |
| `required_files` | `--files` | A path to a file on disk |
| `required_parameters` | `--params` | A string value (e.g., number of threads) |

!!! info "Custom asset classes and recipes"
    To add your own asset types or recipes, see [Data Channels](data_channels.md) and [Flexible Asset Types](flexible_asset_types.md).

## Section 4: Building a single asset

### Build the fasta asset

The fasta asset exports sequences from the RefgetStore to standard FASTA files on disk:

```bash
refgenie1 build rCRSd/fasta
```

Refgenie runs the fasta recipe, which:

1. Exports sequences from the RefgetStore
2. Computes a FASTA index natively (no samtools needed)
3. Writes chromosome sizes
4. Records the asset in the database

Verify the asset is present:

```bash
refgenie1 list
```

Retrieve the path to the asset and its seek keys:

```bash
refgenie1 seek rCRSd/fasta            # default seek key (FASTA file)
refgenie1 seek rCRSd/fasta.fai        # FASTA index
refgenie1 seek rCRSd/fasta.chrom_sizes  # chromosome sizes
```

### Build a derived asset

Build a bowtie2 index. The bowtie2_index recipe requires a fasta asset as input. Refgenie resolves this dependency automatically:

```bash
refgenie1 build rCRSd/bowtie2_index
```

To override the parent asset explicitly:

```bash
refgenie1 build rCRSd/bowtie2_index --assets fasta=rCRSd/fasta:default
```

To run the build inside a Docker container (requires Docker):

```bash
refgenie1 build rCRSd/bowtie2_index -d
```

### Build + stage in one step

The `--stage` flag builds and stages in a single command:

```bash
refgenie1 build rCRSd/bwa_index --stage
```

## Section 5: Staging assets for serving

Staging makes a built asset servable. It creates a `StagedAsset` record in the database and prepares the asset in `genome_stage_folder`.

### Serving modes

Refgenie supports three serving modes, set per asset class. See [Serving modes](serving_modes.md) for the full reference.

| Mode | What is stored in stage folder | Use case |
|------|-------------------------------|---------|
| `file` | Directory symlink to asset in `genome_folder` | Individual file access; files served at their own URLs |
| `archive` | Compressed `.tgz` tarball | Aligner indexes and large bundles; bulk download |
| `none` | Nothing | Metadata-only; users must build locally |

### Stage a file-mode asset

The fasta asset class defaults to `file` mode. Staging creates a symlink in the stage folder:

```bash
refgenie1 stage rCRSd/fasta
```

Inspect the stage folder:

```bash
ls -la /tmp/rg_tutorial/staged/
```

You will see a directory symlink pointing to the built asset in `genome_folder`. Individual files in that directory are each accessible at their own URL when served.

### Stage an archive-mode asset

The bowtie2_index asset class uses `archive` mode. Staging creates a `.tgz` tarball:

```bash
refgenie1 stage rCRSd/bowtie2_index
```

Inspect the stage folder:

```bash
ls /tmp/rg_tutorial/staged/
```

You will see a `.tgz` file. Clients download this tarball and extract it locally.

### List staged assets

```bash
refgenie1 stage staged
```

### Unstage an asset

Unstaging removes the symlink or tarball from the stage folder. The original asset in `genome_folder` is not affected:

```bash
refgenie1 stage unstage rCRSd/fasta
```

## Section 6: Serving assets locally

Start the refgenie server:

```bash
refgenie1 serve
```

By default, this serves on `http://localhost:8000`. To use a different port:

```bash
refgenie1 serve --port 9000
```

### Server endpoints

| Asset mode | Endpoint | Description |
|-----------|---------|-------------|
| Any | `GET /v4/assets` | List all staged assets |
| `file` | `GET /v4/assets/{digest}/files` | List individual files in a file-mode asset |
| `file` | `GET /v4/assets/{digest}/files/{path}` | Download a specific file |
| `archive` | `GET /v4/assets/{digest}/archive` | Download the `.tgz` tarball |
| Any (DRS) | `GET /ga4gh/drs/v1/objects/{id}` | DRS object metadata |

File-mode assets expose each file at its own URL, so downstream tools and web applications can fetch exactly the file they need. Archive-mode assets are served as a single tarball download.

### Remote seek

Get remote URLs to individual files without downloading them:

```bash
refgenie1 seekr rCRSd/fasta -s http://localhost:8000
```

This returns URLs that point directly to the files, suitable for streaming workflows.

### Pull from a client

On another machine that has subscribed to the server, pull assets normally:

```bash
refgenie1 subscribe http://yourserver.example.com
refgenie1 pull rCRSd/fasta
```

Pull behavior adapts to serving mode automatically:

- **archive mode**: downloads the `.tgz` and extracts it
- **file mode**: downloads files individually
- **none mode**: returns an error with instructions to build locally

## Section 7: Building in bulk with PEP and Snakemake

For building many genomes and asset types together, refgenie integrates with [PEP](https://pep.databio.org/) and [Snakemake](https://snakemake.readthedocs.io/).

The workflow:

1. Define what to build in a PEP (sample and subsample tables)
2. Generate a Snakefile from the PEP using refgenie
3. Run snakemake to build all assets in parallel, with correct dependency ordering

### Create a PEP

Create a project directory:

```bash
mkdir -p /tmp/rg_pep
```

**`/tmp/rg_pep/config.yaml`:**

```yaml
pep_version: 2.0.0
sample_table: sample_table.csv
subsample_table: subsample_table.csv

sample_modifiers:
  derive:
    attributes: [fasta_file_path]
    sources:
      local: /tmp/{genome_name}.fa
```

**`/tmp/rg_pep/sample_table.csv`:**

```csv
sample_name,genome_name,genome_description,species_name,fasta_file_path
rCRSd,rCRSd,Revised Cambridge Reference Sequence,Homo sapiens,local
```

**`/tmp/rg_pep/subsample_table.csv`:**

```csv
genome_name,asset_group
rCRSd,fasta
rCRSd,bowtie2_index
rCRSd,bwa_index
```

Each row in `subsample_table.csv` declares one asset to build for a genome. The `asset_group` value must match a recipe name known to refgenie.

### Generate the Snakefile

```bash
refgenie1 generate snakefile --output-path /tmp/rg_pep/build.smk
```

The generated Snakefile contains:

- A `genome_init` rule that initializes each genome from its FASTA file
- Build rules for each asset, derived from the configured recipes
- Stage rules that run after each build
- Dependency chains ensuring that derived assets (e.g., bowtie2_index) are built after their parents (e.g., fasta)

!!! note "Snakefile uses current configuration"
    The Snakefile is generated from your current refgenie configuration — specifically the asset classes and recipes currently registered. If you add new recipes, regenerate the Snakefile.

### Run the workflow

```bash
snakemake --snakefile /tmp/rg_pep/build.smk --jobs 4
```

Snakemake parallelizes across independent assets while respecting dependencies. The `--jobs` flag controls parallel job count.

For cluster execution with SLURM:

```bash
snakemake \
  --snakefile /tmp/rg_pep/build.smk \
  --jobs unlimited \
  --default-resources slurm_account=myaccount slurm_partition=standard mem_mb=32000 \
  --cores 8 \
  --workflow-profile /path/to/snakemake_profile_dir
```

where the profile directory contains a `config.yaml` with `executor: slurm`.

### Verify results

```bash
refgenie1 list
refgenie1 stage staged
```

### Managing software dependencies for recipes

Recipes require bioinformatics tools (e.g., bowtie2-build). Options:

- **Docker** (`refgenie1 build -d`): uses the Docker image specified in the recipe
- **Bulker**: see the [Bulker refgenie tutorial](https://bulker.databio.org/en/latest/refgenie_tutorial/)
- **Snakemake containers**: `snakemake --software-deployment-method apptainer`
- **Native install**: install tools in your PATH

## Section 8: Pushing to cloud remotes

Pushing uploads staged assets to cloud storage (S3, Azure, GCS) so they can be served publicly. Push intent is tracked in the database as `RemoteAssetLink` records, giving you a clear record of what has been pushed and what is pending.

### Configure a remote

```bash
refgenie1 remote add \
  --type s3 \
  --prefix s3://my-bucket/refgenie \
  --name production \
  --push-command "aws s3 cp {local_path} s3://my-bucket/{relative_path}"
```

The `--push-command` template uses these placeholders:

| Placeholder | Value |
|------------|-------|
| `{local_path}` | Absolute path to the staged file or directory |
| `{relative_path}` | Path relative to `genome_stage_folder` |
| `{prefix}` | The remote's prefix string |
| `{genome_stage_folder}` | The `genome_stage_folder` path |

List configured remotes:

```bash
refgenie1 remote list
```

### Push staged assets

Push all staged assets that have not yet been pushed:

```bash
refgenie1 push
```

Preview what would be pushed without uploading:

```bash
refgenie1 push --dry-run
```

Push to a specific remote only:

```bash
refgenie1 push --remote production
```

Push assets for a specific genome only:

```bash
refgenie1 push --genome rCRSd
```

### Check push status

```bash
refgenie1 remote status
```

This shows pushed/unpushed counts per remote along with the digests of unpushed assets.

### Folder sync alternative

For bulk uploads you can use your cloud CLI directly and then mark assets as pushed:

```bash
# Sync the entire stage folder to S3, following symlinks
aws s3 sync /tmp/rg_tutorial/staged/ s3://my-bucket/refgenie/ --follow-symlinks

# Then push using folder_sync strategy to record everything as pushed
refgenie1 push --strategy folder_sync
```

The `folder_sync` strategy marks all staged assets as pushed without re-uploading files, which is useful after a manual sync.

## Section 9: Complete workflow summary

Here is the full lifecycle in a single command sequence:

```bash
# 1. Initialize
refgenie1 init \
  --genome-folder /path/to/genomes \
  --genome-stage-folder /path/to/staged

# 2. Initialize a genome
refgenie1 genome init \
  --fasta rCRSd.fa \
  --name rCRSd \
  --description "Revised Cambridge Reference Sequence"

# 3. Build + stage individual assets
refgenie1 build rCRSd/fasta --stage
refgenie1 build rCRSd/bowtie2_index --stage

# 4. Or build + stage everything in bulk
refgenie1 generate snakefile --output-path build.smk
snakemake --snakefile build.smk --jobs 4

# 5. Serve locally
refgenie1 serve

# 6. Push to cloud
refgenie1 remote add \
  --type s3 \
  --prefix s3://bucket/refgenie \
  --name production \
  --push-command "aws s3 cp {local_path} s3://bucket/{relative_path}"
refgenie1 push
```

### Quick reference

| Command | Description |
|---------|-------------|
| `refgenie1 init` | Initialize refgenie configuration and database |
| `refgenie1 genome init --fasta FILE --name NAME` | Register a genome from a FASTA file |
| `refgenie1 genome list` | List all registered genomes |
| `refgenie1 asset_class list` | List available asset classes |
| `refgenie1 recipe list` | List available recipes |
| `refgenie1 recipe show NAME -q` | Show recipe requirements |
| `refgenie1 build GENOME/ASSET` | Build an asset |
| `refgenie1 build GENOME/ASSET --stage` | Build and stage in one step |
| `refgenie1 list` | List all local assets |
| `refgenie1 seek GENOME/ASSET` | Get local file path |
| `refgenie1 stage GENOME/ASSET` | Stage a built asset |
| `refgenie1 stage unstage GENOME/ASSET` | Unstage an asset |
| `refgenie1 stage staged` | List all staged assets |
| `refgenie1 serve` | Start the local refgenie server |
| `refgenie1 seekr GENOME/ASSET -s URL` | Get remote file URL |
| `refgenie1 generate snakefile -o FILE` | Generate a Snakemake build workflow |
| `refgenie1 remote add` | Configure a cloud push remote |
| `refgenie1 remote list` | List configured remotes |
| `refgenie1 remote status` | Show push status per remote |
| `refgenie1 push` | Push staged assets to cloud remotes |
| `refgenie1 push --dry-run` | Preview push without uploading |
