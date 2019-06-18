# Usage reference

## `refgenie --help`

```console
version: 0.4.1
usage: refgenie [-h] [-V] [--logdev] [--verbosity V] [--silent]
                {pull,seek,list,build,init,listr} ...

refgenie - builds and manages reference genome assemblies

positional arguments:
  {pull,seek,list,build,init,listr}
    pull                Download assets.
    seek                Get the path to a local asset
    list                List available local genomes.
    build               Build genome assets
    init                Initialize a genome configuration.
    listr               List available genomes and assets on server.

optional arguments:
  -h, --help            show this help message and exit
  -V, --version         show program's version number and exit
  --logdev              Expand content of logging message format.
  --verbosity V         Set logging level (1-5 or logging module level name)
  --silent              Silence logging. Overrides --verbosity.

https://refgenie.databio.org
```

## `refgenie init --help`

```console
version: 0.4.1
usage: refgenie init [-h] -c GENOME_CONFIG [-s GENOME_SERVER]

Initialize a genome configuration.

optional arguments:
  -h, --help            show this help message and exit
  -c GENOME_CONFIG, --genome-config GENOME_CONFIG
                        Path to local genome configuration file, to read from
                        and/or to create or update, depending on the operation
  -s GENOME_SERVER, --genome-server GENOME_SERVER
                        URL to use for the genome_server attribute in config
                        file. Defaults : http://refgenomes.databio.org
```

## `refgenie list --help`

```console
version: 0.4.1
usage: refgenie list [-h] [-c GENOME_CONFIG]

List available local genomes.

optional arguments:
  -h, --help            show this help message and exit
  -c GENOME_CONFIG, --genome-config GENOME_CONFIG
                        Path to local genome configuration file, to read from
                        and/or to create or update, depending on the operation
```

## `refgenie listr --help`

```console
version: 0.4.1
usage: refgenie listr [-h] [-c GENOME_CONFIG]

List available genomes and assets on server.

optional arguments:
  -h, --help            show this help message and exit
  -c GENOME_CONFIG, --genome-config GENOME_CONFIG
                        Path to local genome configuration file, to read from
                        and/or to create or update, depending on the operation
```

## `refgenie pull --help`

```console
version: 0.4.1
usage: refgenie pull [-h] [-c GENOME_CONFIG] -g GENOME -a ASSET [ASSET ...]
                     [-u]

Download assets.

optional arguments:
  -h, --help            show this help message and exit
  -c GENOME_CONFIG, --genome-config GENOME_CONFIG
                        Path to local genome configuration file, to read from
                        and/or to create or update, depending on the operation
  -g GENOME, --genome GENOME
                        Reference assembly ID, e.g. mm10
  -a ASSET [ASSET ...], --asset ASSET [ASSET ...]
                        Name of asset, a key in a genome config file
  -u, --no-untar        Do not extract tarballs.
```

