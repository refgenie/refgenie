# Usage reference

## `refgenie --help`

```console
version: 0.4.0-dev
usage: refgenie [-h] [-V] [--silent] [--logdev] [--verbosity V]
                {pull,init,build,seek,list,listr} ...

refgenie - builds and manages reference genome assemblies

positional arguments:
  {pull,init,build,seek,list,listr}
    pull                Download assets.
    init                Initialize a genome configuration.
    build               Build genome assets
    seek                Get the path to a local asset
    list                List available local genomes.
    listr               List available genomes and assets on server.

optional arguments:
  -h, --help            show this help message and exit
  -V, --version         show program's version number and exit
  --silent              Silence logging. Overrides --verbosity.
  --logdev              Expand content of logging message format.
  --verbosity V         Set logging level (1-5 or logging module level name)

https://refgenie.databio.org
```

