# How to start with refgenie (by Alex)

### 1. Install refgenie

```python
uv pip install . # or pip install refgenie # python should be >=3.11
```

### 2. Create db config files. eg.

```yaml
type: sqlite
name: /home/bnt4me/virginia/refgenie/config/refgenie.sqlite
```

** Note:** If you don't have a config file, refgenie will create one for you in the `.refgenie` directory.
** Note2:** If you provided incorrect path, refgenie will create a new file in the `.refgenie` directory/ or throw an error if the directory does not exist.

### 3. Set all necessary env vars:

```bash
export REFGENIE_GENOME_STAGE_FOLDER=/home/bnt4me/virginia/refgenie/archive/
export REFGENIE_GENOME_FOLDER=/home/bnt4me/virginia/refgenie/genomes
export REFGENIE_DB_CONFIG_PATH=/home/bnt4me/virginia/refgenie/config/sqlite_db.yaml
export REFGENIE_INPUTS=/home/bnt4me/virginia/refgenie/inputs
```

### 4. Initialize refgenie

```bash
refgenie init
```

This function will create all required folders and pull default asset class: `fasta`.


## 5. Connect to the data channel.

We need one data channel where we can fetch all asset classes, to be able later to build assets.

```python
refgenie data_channel add my-fav-channel https https://refgenie.github.io/refgenie-registry/index.yaml
```

### 6. List your data channels

```bash
refgenie data_channel list
```

```text
                                              Data Channels                                               
┏━━━━━━━━━━━━━━━━┳━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━┓
┃ Name           ┃ Type  ┃ Index Address                                 ┃ Description ┃ Credentials set ┃
┡━━━━━━━━━━━━━━━━╇━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━┩
│ my-fav-channel │ https │ https://refgenie.github.io/refgenie-registry/index.yaml │             │ False           │
└────────────────┴───────┴───────────────────────────────────────────────┴─────────────┴─────────────────┘
```

### 7. Sync asset classes and recipes from the data channel

```bash
refgenie data_channel sync my-fav-channel --exists-ok
```

Now all the asset classes and recipes are registered in the refgenie database.

Another way to add asset classes and recipes is to use local path:

```bash

refgenie asset_class add --source path/to/asset_class.yaml
refgenie recipe add --source path/to/recipe.yaml

```


### 8. List all asset classes:

```bash
refgenie asset_class list
```

Example output:

```text
                                                                         Asset Classes                                                                           
┏━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Name                    ┃ Version ┃ Seek keys                                                   ┃ Description                                                  ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ fasta                   │ 0.0.1   │ chrom_sizes, fai, fasta                                     │ Sequences in the FASTA format, indexed FASTA (produced with  │
│                         │         │                                                             │ samtools index) and chromosome sizes file                    │
├─────────────────────────┼─────────┼─────────────────────────────────────────────────────────────┼──────────────────────────────────────────────────────────────┤
│ abundant_sequences      │ 0.0.1   │ abundant_sequences, adapter_contam, phix, polyA, polyC      │ Abundant sequences in the FASTA format -- PhiX spike-in,     │
│                         │         │                                                             │ Poly(A), Poly(C) and adapter sequences                       │
├─────────────────────────┼─────────┼─────────────────────────────────────────────────────────────┼──────────────────────────────────────────────────────────────┤
│ bed                     │ 0.0.1   │ bedgz                                                       │ Genomic feature annotations asset which provides access to   │
│                         │         │                                                             │ all annotated transcripts in BED format                      │
├─────────────────────────┼─────────┼─────────────────────────────────────────────────────────────┼──────────────────────────────────────────────────────────────┤
│
```

### 8. List all asset classes:

```bash
refgenie recipe list
```

Example output:

```text
┏━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Name                    ┃ Version ┃ Seek keys                                                                 ┃ Description                                                                                                      ┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ fasta                   │ 0.0.1   │ chrom_sizes, fai, fasta                                                   │ Sequences in the FASTA format, indexed FASTA (produced with samtools index) and chromosome sizes file            │
├─────────────────────────┼─────────┼───────────────────────────────────────────────────────────────────────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ abundant_sequences      │ 0.0.1   │ abundant_sequences, adapter_contam, phix, polyA, polyC                    │ Abundant sequences in the FASTA format -- PhiX spike-in, Poly(A), Poly(C) and adapter sequences                  │
├─────────────────────────┼─────────┼───────────────────────────────────────────────────────────────────────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
│ bed                     │ 0.0.1   │ bedgz                                                                     │ Genomic feature annotations asset which provides access to all annotated transcripts in BED format               │
├─────────────────────────┼─────────┼───────────────────────────────────────────────────────────────────────────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤
```



# How to build assets:
To build asset you need to have all necessary files + required recopies and asset classes registered in the refgenie database.

e.g. if we have `baw_index` asset class and recipe, we can build a BWA index asset with the following command:

```bash
refgenie build rCRSd/bwa_index
```


# How to add fasta files as asset for new organism: 

```bash

refgenie build bBucAby1.pri/fasta:default --genome-description "Abyssinian ground-hornbill (primary hap 2019)" --files fasta=/home/bnt4me/virginia/refgenie/my_fastas/GCA_009769605.1.fa.gz --archive

```

# How to use bulker for all the requirements

1. Initialize bulker locally:

Set bash variables and init bulker:
```bash
rm "bulker_config.yaml"
export BULKERCFG="/home/bnt4me/.bulker_config.yaml"
bulker init -c $BULKERCFG # DO it if not initialized yet
```

2. Load your bulker

```bash
bulker load refgenie -m /home/bnt4me/virginia/refgenie/refgenie_bulker_manifest.yaml
```

3. Activate bulker

```bash
eval "$(bulker activate refgenie -e)"
```

# Vocabularies:

- **asset class**: 
- **recipe**: 
- **data channel**: Refgenie data channels provide asset class and recipe definitions.
- **asset**: A specific instance of an asset class built using a recipe.
- **alias**: A human-readable name for an asset, often used to simplify access to the asset.