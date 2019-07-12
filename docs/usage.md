# Usage reference

## `refgenie --help`

```console
version: 0.5.0
usage: refgenie [-h] [-V] [--verbosity V] [--logdev] [--silent]
                {listr,add,list,build,pull,seek,init} ...

refgenie - builds and manages reference genome assemblies

positional arguments:
  {listr,add,list,build,pull,seek,init}
    listr               List available remote assets.
    add                 Add local asset to the config file.
    list                List available local assets.
    build               Build genome assets.
    pull                Download assets.
    seek                Get the path to a local asset.
    init                Initialize a genome configuration.

optional arguments:
  -h, --help            show this help message and exit
  -V, --version         show program's version number and exit
  --verbosity V         Set logging level (1-5 or logging module level name)
  --logdev              Expand content of logging message format.
  --silent              Silence logging. Overrides --verbosity.

https://refgenie.databio.org
```

## `refgenie init --help`

```console
version: 0.5.0
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
version: 0.5.0
usage: refgenie list [-h] [-c GENOME_CONFIG]

List available local assets.

optional arguments:
  -h, --help            show this help message and exit
  -c GENOME_CONFIG, --genome-config GENOME_CONFIG
                        Path to local genome configuration file.
```

## `refgenie seek --help`

```console
version: 0.5.0
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
version: 0.5.0
usage: refgenie listr [-h] [-c GENOME_CONFIG]

List available remote assets.

optional arguments:
  -h, --help            show this help message and exit
  -c GENOME_CONFIG, --genome-config GENOME_CONFIG
                        Path to local genome configuration file.
```

## `refgenie pull --help`

```console
version: 0.5.0
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
version: 0.5.0
usage: refgenie add [-h] [-c GENOME_CONFIG] -g GENOME -a ASSET [ASSET ...] -p
                    PATH

Add local asset to the config file.

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

