# Usage reference

## `refgenie --help`

```console
version: 0.7.3
usage: refgenie [-h] [--version] [--silent] [--verbosity V] [--logdev]
                {init,list,listr,pull,build,seek,add,remove,getseq,tag,id,subscribe,unsubscribe}
                ...

refgenie - reference genome asset manager

positional arguments:
  {init,list,listr,pull,build,seek,add,remove,getseq,tag,id,subscribe,unsubscribe}
    init                Initialize a genome configuration.
    list                List available local assets.
    listr               List available remote assets.
    pull                Download assets.
    build               Build genome assets.
    seek                Get the path to a local asset.
    add                 Add local asset to the config file.
    remove              Remove a local asset.
    getseq              Get sequences from a genome.
    tag                 Tag an asset.
    id                  Return the asset digest.
    subscribe           Add a refgenieserver URL to the config.
    unsubscribe         Remove a refgenieserver URL from the config

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --silent              Silence logging. Overrides verbosity.
  --verbosity V         Set logging level (1-5 or logging module level name)
  --logdev              Expand content of logging message format.

https://refgenie.databio.org

```

## `refgenie init --help`

```console
usage: refgenie init [-h] -c GENOME_CONFIG [-s GENOME_SERVER]

Initialize a genome configuration.

optional arguments:
  -h, --help            show this help message and exit
  -c GENOME_CONFIG, --genome-config GENOME_CONFIG
                        Path to local genome configuration file. Optional if
                        REFGENIE environment variable is set.
  -s GENOME_SERVER, --genome-server GENOME_SERVER
                        URL to use for the genome_servers attribute in config
                        file. Default: http://refgenomes.databio.org
```

## `refgenie list --help`

```console
usage: refgenie list [-h] [-c GENOME_CONFIG] [-g GENOME]

List available local assets.

optional arguments:
  -h, --help            show this help message and exit
  -c GENOME_CONFIG, --genome-config GENOME_CONFIG
                        Path to local genome configuration file. Optional if
                        REFGENIE environment variable is set.
  -g GENOME, --genome GENOME
                        Reference assembly ID, e.g. mm10
```

## `refgenie listr --help`

```console
usage: refgenie listr [-h] [-c GENOME_CONFIG] [-g GENOME]

List available remote assets.

optional arguments:
  -h, --help            show this help message and exit
  -c GENOME_CONFIG, --genome-config GENOME_CONFIG
                        Path to local genome configuration file. Optional if
                        REFGENIE environment variable is set.
  -g GENOME, --genome GENOME
                        Reference assembly ID, e.g. mm10
```

## `refgenie pull --help`

```console
usage: refgenie pull [-h] [-c GENOME_CONFIG] [-g GENOME] [-u]
                     asset-registry-paths [asset-registry-paths ...]

Download assets.

positional arguments:
  asset-registry-paths  One or more registry path strings that identify assets
                        (e.g. hg38/fasta or hg38/fasta:tag)

optional arguments:
  -h, --help            show this help message and exit
  -c GENOME_CONFIG, --genome-config GENOME_CONFIG
                        Path to local genome configuration file. Optional if
                        REFGENIE environment variable is set.
  -g GENOME, --genome GENOME
                        Reference assembly ID, e.g. mm10
  -u, --no-untar        Do not extract tarballs.
```

## `refgenie build --help`

```console
usage: refgenie build [-h] [-c GENOME_CONFIG] [-R] [-C CONFIG_FILE] [-N]
                      [--tag-description TAG_DESCRIPTION]
                      [--genome-description GENOME_DESCRIPTION] [-d]
                      [--assets ASSETS [ASSETS ...]]
                      [--files FILES [FILES ...]]
                      [--params PARAMS [PARAMS ...]]
                      [-v VOLUMES [VOLUMES ...]] [-o OUTFOLDER] [-q]
                      [-r RECIPE] [-g GENOME]
                      asset-registry-paths [asset-registry-paths ...]

Build genome assets.

positional arguments:
  asset-registry-paths  One or more registry path strings that identify assets
                        (e.g. hg38/fasta or hg38/fasta:tag)

optional arguments:
  -h, --help            show this help message and exit
  -c GENOME_CONFIG, --genome-config GENOME_CONFIG
                        Path to local genome configuration file. Optional if
                        REFGENIE environment variable is set.
  -R, --recover         Overwrite locks to recover from previous failed run
  -C CONFIG_FILE, --config CONFIG_FILE
                        Pipeline configuration file (YAML). Relative paths are
                        with respect to the pipeline script.
  -N, --new-start       Overwrite all results to start a fresh run
  --tag-description TAG_DESCRIPTION
                        Add tag level description (e.g. built with version
                        0.3.2)
  --genome-description GENOME_DESCRIPTION
                        Add genome level description (e.g. The mouse
                        mitochondrial genome, released in Dec 2013)
  -d, --docker          Run all commands in the refgenie docker container.
  --assets ASSETS [ASSETS ...]
                        Override the default genome, asset and tag of the
                        parents (e.g. fasta=hg38/fasta:default
                        gtf=mm10/gencode_gtf:default).
  --files FILES [FILES ...]
                        Provide paths to the required files (e.g.
                        fasta=/path/to/file.fa.gz).
  --params PARAMS [PARAMS ...]
                        Provide required parameter values (e.g.
                        param1=value1).
  -v VOLUMES [VOLUMES ...], --volumes VOLUMES [VOLUMES ...]
                        If using docker, also mount these folders as volumes.
  -o OUTFOLDER, --outfolder OUTFOLDER
                        Override the default path to genomes folder, which is
                        the genome_folder attribute in the genome
                        configuration file.
  -q, --requirements    Show the build requirements for the specified asset
                        and exit.
  -r RECIPE, --recipe RECIPE
                        Provide a recipe to use.
  -g GENOME, --genome GENOME
                        Reference assembly ID, e.g. mm10
```

## `refgenie seek --help`

```console
usage: refgenie seek [-h] [-c GENOME_CONFIG] [-g GENOME]
                     asset-registry-paths [asset-registry-paths ...]

Get the path to a local asset.

positional arguments:
  asset-registry-paths  One or more registry path strings that identify assets
                        (e.g. hg38/fasta or hg38/fasta:tag or
                        hg38/fasta.fai:tag)

optional arguments:
  -h, --help            show this help message and exit
  -c GENOME_CONFIG, --genome-config GENOME_CONFIG
                        Path to local genome configuration file. Optional if
                        REFGENIE environment variable is set.
  -g GENOME, --genome GENOME
                        Reference assembly ID, e.g. mm10
```

## `refgenie add --help`

```console
usage: refgenie add [-h] [-c GENOME_CONFIG] [-g GENOME] -p PATH
                    asset-registry-paths [asset-registry-paths ...]

Add local asset to the config file.

positional arguments:
  asset-registry-paths  One or more registry path strings that identify assets
                        (e.g. hg38/fasta or hg38/fasta:tag)

optional arguments:
  -h, --help            show this help message and exit
  -c GENOME_CONFIG, --genome-config GENOME_CONFIG
                        Path to local genome configuration file. Optional if
                        REFGENIE environment variable is set.
  -g GENOME, --genome GENOME
                        Reference assembly ID, e.g. mm10
  -p PATH, --path PATH  Relative local path to asset
```

## `refgenie remove --help`

```console
version: 0.5.0
usage: refgenie remove [-h] [-c GENOME_CONFIG] -g GENOME
                       [-a ASSET [ASSET ...]]

Remove a local asset.

optional arguments:
  -h, --help            show this help message and exit
  -c GENOME_CONFIG, --genome-config GENOME_CONFIG
                        Path to local genome configuration file.
  -g GENOME, --genome GENOME
                        Reference assembly ID, e.g. mm10
  -a ASSET [ASSET ...], --asset ASSET [ASSET ...]
                        Name of one or more assets (keys in genome config
                        file)
```

## `refgenie getseq --help`

```console
usage: refgenie getseq [-h] [-c GENOME_CONFIG] -g GENOME -l LOCUS

Get sequences from a genome.

optional arguments:
  -h, --help            show this help message and exit
  -c GENOME_CONFIG, --genome-config GENOME_CONFIG
                        Path to local genome configuration file. Optional if
                        REFGENIE environment variable is set.
  -g GENOME, --genome GENOME
                        Reference assembly ID, e.g. mm10
  -l LOCUS, --locus LOCUS
                        Coordinates to retrieve sequence for; such has
                        'chr1:50000-50200'.
```

## `refgenie tag --help`

```console
usage: refgenie tag [-h] [-c GENOME_CONFIG] [-g GENOME] -t TAG
                    asset-registry-paths [asset-registry-paths ...]

Assign a selected tag to an asset.

positional arguments:
  asset-registry-paths  One or more registry path strings that identify assets
                        (e.g. hg38/fasta or hg38/fasta:tag)

optional arguments:
  -h, --help            show this help message and exit
  -c GENOME_CONFIG, --genome-config GENOME_CONFIG
                        Path to local genome configuration file. Optional if
                        REFGENIE environment variable is set.
  -g GENOME, --genome GENOME
                        Reference assembly ID, e.g. mm10
  -t TAG, --tag TAG     Tag to assign to an asset
```


## `refgenie id --help`

```console
usage: refgenie id [-h] [-c GENOME_CONFIG] [-g GENOME]
                   asset-registry-paths [asset-registry-paths ...]

Return the asset digest.

positional arguments:
  asset-registry-paths  One or more registry path strings that identify assets
                        (e.g. hg38/fasta or hg38/fasta:tag)

optional arguments:
  -h, --help            show this help message and exit
  -c GENOME_CONFIG, --genome-config GENOME_CONFIG
                        Path to local genome configuration file. Optional if
                        REFGENIE environment variable is set.
  -g GENOME, --genome GENOME
                        Reference assembly ID, e.g. mm10
```

## `refgenie subscribe --help`

```console
usage: refgenie subscribe [-h] [-c GENOME_CONFIG] [-r] -s GENOME_SERVER
                          [GENOME_SERVER ...]

Add a refgenieserver URL to the config.

optional arguments:
  -h, --help            show this help message and exit
  -c GENOME_CONFIG, --genome-config GENOME_CONFIG
                        Path to local genome configuration file. Optional if
                        REFGENIE environment variable is set.
  -r, --reset           Overwrite the current list of server URLs
  -s GENOME_SERVER [GENOME_SERVER ...], --genome-server GENOME_SERVER [GENOME_SERVER ...]
                        One or more URLs to add to the genome_servers
                        attribute in config file
```


## `refgenie unsubscribe --help`

```console
usage: refgenie unsubscribe [-h] [-c GENOME_CONFIG] -s GENOME_SERVER
                            [GENOME_SERVER ...]

Remove a refgenieserver URL from the config

optional arguments:
  -h, --help            show this help message and exit
  -c GENOME_CONFIG, --genome-config GENOME_CONFIG
                        Path to local genome configuration file. Optional if
                        REFGENIE environment variable is set.
  -s GENOME_SERVER [GENOME_SERVER ...], --genome-server GENOME_SERVER [GENOME_SERVER ...]
                        One or more URLs to add to the genome_servers
                        attribute in config file
```