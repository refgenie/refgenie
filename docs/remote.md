# Remote mode in refgenie

Starting with version 0.11.0, refgenie can be used in *remote mode*, which means that in some cases the genome configuration file is not required. **Therefore, you can skip `refgenie init` and start workling with refgenie right after installation!**

## Commands available in remote mode

*Hint: all of these commands end with "r"*

There are a few commands that do not require genome configuration file to run. `-s/--genome-servers` argument specifies the list of servers you want refgenie to query. Default server ([http://refgenomes.databio.org](http://refgenomes.databio.org)) is used if not provided.

### List remote assets with `refgenie listr`

You can list assets available on remote servers with `refgenie listr`.

```console
~ refgenie listr -s http://rg.databio.org
Using default config. No config found in env var: ['REFGENIE']
Subscribed to: http://rg.databio.org

                        Remote refgenie assets
                  Server URL: http://rg.databio.org
┏━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ genome              ┃ assets                                       ┃
┡━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ mouse_chrM2x        │ fasta, bowtie2_index, bwa_index              │
│ rCRSd               │ fasta, bowtie2_index                         │
│ human_repeats       │ fasta, hisat2_index, bwa_index               │
│ hg38                │ fasta, bowtie2_index                         │
└─────────────────────┴──────────────────────────────────────────────┘
        use refgenie listr -g <genome> for more detailed view
```

### Find remote asset paths with `refgenie seekr`

You can seek for remote asset paths with `refgenie seekr`:

```console
~ refgenie seekr hg38/fasta -s http://rg.databio.org -r s3

Using default config. No config found in env var: ['REFGENIE']
Subscribed to: http://rg.databio.org
No local digest for genome alias: hg38
Setting 'hg38' identity with server: http://rg.databio.org/v3/genomes/genome_digest/hg38
Determined server digest for local genome alias (hg38): 2230c535660fb4774114bfa966a62f823fdb6d21acf138d4
Set genome alias (2230c535660fb4774114bfa966a62f823fdb6d21acf138d4: hg38)
s3://awspds.refgenie.databio.org/rg.databio.org/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/fasta__default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.fa
```

`-r`/`--remote-class` command specifies the data provider link to be used in the output. Please refer to the refgenieserver instance API endpoint that lists available options, for example: [http://rg.databio.org/remotes/dict](http://rg.databio.org/remotes/dict).

### Replace asset registry paths with remote asset paths using `refgenie populater`

You can replace refgenie asset registry paths in **text** or **files** with `refgenie populater`. Any string that matches the following format will be replaced with a remote path:

```console
refgenie://genome_alias/asset.seek_key:tag
```


#### populate text from standard input

```console
~ echo 'test remote populating refgenie://hg38/fasta.fasta:default' | refgenie populater -s http://rg.databio.org -r s3

Using default config. No config found in env var: ['REFGENIE']
Subscribed to: http://rg.databio.org
No local digest for genome alias: hg38
Setting 'hg38' identity with server: http://rg.databio.org/v3/genomes/genome_digest/hg38
Determined server digest for local genome alias (hg38): 2230c535660fb4774114bfa966a62f823fdb6d21acf138d4
Set genome alias (2230c535660fb4774114bfa966a62f823fdb6d21acf138d4: hg38)
test remote populating s3://awspds.refgenie.databio.org/rg.databio.org/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/fasta__default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.fa:default
```

`-r`/`--remote-class` command specifies the data provider link to be used in the output. Please refer to the refgenieserver instance API endpoint that lists available options, for example: [http://rg.databio.org/remotes/dict](http://rg.databio.org/remotes/dict).

#### populate text from file

- Check input file contents
```console
~ cat remote_populate_test.txt

human genome FASTA file: refgenie://hg38/fasta.fasta
yeast doubled genome FASTA file: refgenie://rCRSd/fasta.fasta
```

- Run `refgenie populater`
```console
~ refgenie populater -f remote_populate_test.txt -s http://rg.databio.org > remote_populate_test_output.txt

Using default config. No config found in env var: ['REFGENIE']
Subscribed to: http://rg.databio.org
No local digest for genome alias: hg38
Setting 'hg38' identity with server: http://rg.databio.org/v3/genomes/genome_digest/hg38
Determined server digest for local genome alias (hg38): 2230c535660fb4774114bfa966a62f823fdb6d21acf138d4
Set genome alias (2230c535660fb4774114bfa966a62f823fdb6d21acf138d4: hg38)
No local digest for genome alias: rCRSd
Setting 'rCRSd' identity with server: http://rg.databio.org/v3/genomes/genome_digest/rCRSd
Determined server digest for local genome alias (rCRSd): 94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4
Set genome alias (94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4: rCRSd)
```

- Check output file contents
```console
~ cat remote_populate_test_output.txt

human genome FASTA file: http://awspds.refgenie.databio.org/rg.databio.org/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/fasta__default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.fa
yeast doubled genome FASTA file: http://awspds.refgenie.databio.org/rg.databio.org/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/fasta__default/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4.fa
```

## Motivation

The motivation behind the remote mode in refgenie is *cloud computing*. It is becoming a common practice to farm out jobs that require refgenie assets to computing clusters, where refgenie environment is not configured.

Up until now, the user was expected to `init` the refgenie config, `pull` desired assets and then `seek` the path in order to pass it to the data processing workflow. With the new `seekr` command the configuration and data acquisition steps can be skipped. What is more, refgenieserver software provides full flexibility regarding place where the asset files and archives are stored. Therefore, in some cases the data may be readily available within the cloud services provider's servers. For example, [http://refgenomes.databio.org](http://refgenomes.databio.org) refgenieserver instance stores the data in AWS S3, so any jobs running on AWS servers would benefit from the increased performance.
