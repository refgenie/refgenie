# Remote mode in refgenie for cloud computing use

Starting with version 0.11.0, refgenie offers a series of new functions that work directly with remote assets, which we refer to as *remote mode*. Remote mode simply means that we omit the local configuration file and storage of assets, instead relying only on the server for both metadata and asset files. **With remote mode, you can skip `refgenie init` and start workling with refgenie right after installation!**

## Rationale

The motivation behind the remote mode in refgenie is cloud computing. It is becoming a common practice to farm out jobs that require refgenie assets to cloud instances. Remote mode offers 3 key features to facilitate this:

### 1. No initialization required.

Up until now, the user was expected to `init` the refgenie config, `pull` desired assets and then `seek` the path in order to pass it to the data processing workflow. The local config would keep track of which assets had been pulled. This is still the best way to use refgenie when dealing with local computing environments, like a desktop or on-premesis HPC compute cluster. But for cloud computing, with the new `seekr` command, the `init` step may be skipped, and you can simply rely on the cloud server's configuration without any local configuration file needed.

### 2. Unarchived assets are now available on the server

For local computing, users will `pull` assets, which downloads and unarchives them for local computing. For cloud computing, the `pull` step can be skipped; instead, you can now refer to the *unarchived* data directly on the cloud server. For example, [http://refgenomes.databio.org](http://refgenomes.databio.org) refgenieserver instance stores the data in AWS S3, so any jobs running on AWS servers can simply use the unarchived files available direclty on S3, referring to individual files rather than complete archives.

### 3. Flexibilty to choose either HTTP or S3 links.

Refgenieserver instances may now specify multiple protocols, which can be requested by the user. For example, our default server provides either an `http://` URI or an `s3://` URI. This way, if your tooling requires an s3-compatible interface, you can also use refgenie easily by specifying that you want `s3`-compatible URIs. Developers who run their own refgenieserver instances have complete control over the number and type of URIs that can be queried by users.

## Commands available in remote mode

*Hint: all remote mode commands end with "r"*

The new remote commands do not require genome configuration file to run. Instead, the `-s/--genome-servers` argument specifies the list of servers you want refgenie to query. The default server ([http://refgenomes.databio.org](http://refgenomes.databio.org)) is used if `-s` is not provided. The remote mode commands are: `listr`; `seekr`; and `populater`.

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

Since refgenie server stores unarchived assets directly for cloud computing, you can seek for direct remote asset paths with `refgenie seekr`:

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

The path returned to `stdout` is the final line above, `s3://...`, which specifies a direct path to an uncompressed fasta file on the server. Notice we used
`-r`/`--remote-class`, which here specifies the data provider link to be used in the output. Please refer to the refgenieserver instance API endpoint that lists available options, for example: [http://rg.databio.org/remotes/dict](http://rg.databio.org/remotes/dict). For example. you could use `-r http` to get an HTTP link to the file, or `-r s3` to get an S3 URI.

### Replace asset registry paths with remote asset paths using `refgenie populater`

One useful recent feature is to replace refgenie asset registry paths in **text** or **files** with `populate`, which is described in [guide to refgenie populate](populate.md). In remote mode, you use `refgenie populater` (think of it as `refgenie populate` Remote). Just like `populate`, `populatr` will convert any string that matches the refgenie registry path format with a remote path. For example:

```console
refgenie://genome_alias/asset.seek_key:tag
```

You can replace this with `s3://...`. Here are a few real examples:

#### Example 1: populate text from standard input

Here we'll echo some text to stdin, and then use `-r s3` with `refgenie populater` to get the remote s3 link to the specified asset:

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

#### Example 2: populate text in file

We can also use `-f` to specify an input file, which may contain refgenie registry paths. Here's an example file that includes some refgenie paths:

```console
~ cat remote_populate_test.txt

human genome FASTA file: refgenie://hg38/fasta.fasta
doubled mtDNA FASTA file: refgenie://rCRSd/fasta.fasta
```

Run `refgenie populater` to convert those paths into http paths (here, we're not passing `-r s3` so it will return the default links, which is `http`).

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

Here are the resulting file contents:

```console
~ cat remote_populate_test_output.txt

human genome FASTA file: http://awspds.refgenie.databio.org/rg.databio.org/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/fasta__default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.fa
doubled mtDNA FASTA file: http://awspds.refgenie.databio.org/rg.databio.org/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4/fasta__default/94e0d21feb576e6af61cd2a798ad30682ef2428bb7eabbb4.fa
```
