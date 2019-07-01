# Usage reference

## `refgenie --help`

```console
version: 0.4.5-dev
usage: refgenie [-h] [-V] [--silent] [--logdev] [--verbosity V]
                {pull,add,build,listr,list,init,seek} ...

refgenie - builds and manages reference genome assemblies

positional arguments:
  {pull,add,build,listr,list,init,seek}
    pull                Download assets.
    add                 Insert a local asset into the configuration file.
    build               Build genome assets.
    listr               List available genomes and assets on server.
    list                List available local genomes.
    init                Initialize a genome configuration.
    seek                Get the path to a local asset.

optional arguments:
  -h, --help            show this help message and exit
  -V, --version         show program's version number and exit
  --silent              Silence logging. Overrides --verbosity.
  --logdev              Expand content of logging message format.
  --verbosity V         Set logging level (1-5 or logging module level name)

https://refgenie.databio.org
```

## `refgenie init --help`

```console
version: 0.4.5-dev
usage: refgenie init [-h] -c GENOME_CONFIG [-s GENOME_SERVER]

Initialize a genome configuration.

optional arguments:
  -h, --help            show this help message and exit
  -c GENOME_CONFIG, --genome-config GENOME_CONFIG
                        Path to local genome configuration file.
  -s GENOME_SERVER, --genome-server GENOME_SERVER
                        URL to use for the genome_server attribute in config
                        file. Defaults : http://refgenomes.databio.org
```

## `refgenie list --help`

```console
version: 0.4.5-dev
usage: refgenie list [-h] [-c GENOME_CONFIG]

List available local genomes.

optional arguments:
  -h, --help            show this help message and exit
  -c GENOME_CONFIG, --genome-config GENOME_CONFIG
                        Path to local genome configuration file.
```

## `refgenie seek --help`

```console
version: 0.4.5-dev
usage: refgenie seek [-h] [-c GENOME_CONFIG] -g GENOME -a ASSET [ASSET ...]

Get the path to a local asset.

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

## `refgenie listr --help`

```console
version: 0.4.5-dev
usage: refgenie listr [-h] [-c GENOME_CONFIG]

List available genomes and assets on server.

optional arguments:
  -h, --help            show this help message and exit
  -c GENOME_CONFIG, --genome-config GENOME_CONFIG
                        Path to local genome configuration file.
```

## `refgenie pull --help`

```console
version: 0.4.5-dev
usage: refgenie pull [-h] [-c GENOME_CONFIG] -g GENOME -a ASSET [ASSET ...]
                     [-u]

Download assets.

optional arguments:
  -h, --help            show this help message and exit
  -c GENOME_CONFIG, --genome-config GENOME_CONFIG
                        Path to local genome configuration file.
  -g GENOME, --genome GENOME
                        Reference assembly ID, e.g. mm10
  -a ASSET [ASSET ...], --asset ASSET [ASSET ...]
                        Name of one or more assets (keys in genome config
                        file)
  -u, --no-untar        Do not extract tarballs.
```

## `refgenie add --help`

```console
version: 0.4.5-dev
usage: refgenie add [-h] [-c GENOME_CONFIG] -g GENOME -a ASSET [ASSET ...] -p
                    PATH

Insert a local asset into the configuration file.

optional arguments:
  -h, --help            show this help message and exit
  -c GENOME_CONFIG, --genome-config GENOME_CONFIG
                        Path to local genome configuration file.
  -g GENOME, --genome GENOME
                        Reference assembly ID, e.g. mm10
  -a ASSET [ASSET ...], --asset ASSET [ASSET ...]
                        Name of one or more assets (keys in genome config
                        file)
  -p PATH, --path PATH  Relative path to asset
```

