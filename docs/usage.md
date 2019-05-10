# Usage reference

## `refgenie --help`

```console
usage: refgenie [-h] [-R] [-C CONFIG_FILE] -i INPUT [-n NAME] [-a ANNOTATION]
                [-d] [-o OUTFOLDER]

Pipeline

optional arguments:
  -h, --help            show this help message and exit
  -R, --recover         Overwrite locks to recover from previous failed run
  -C CONFIG_FILE, --config CONFIG_FILE
                        Pipeline configuration file (YAML). Relative paths are
                        with respect to the pipeline script.
  -i INPUT, --input INPUT
                        Local path or URL to genome sequence file in .fa,
                        .fa.gz, or .2bit format.
  -n NAME, --name NAME  Name of the genome to build. If ommitted, refgenie
                        will usethe basename of the file specified in --input
  -a ANNOTATION, --annotation ANNOTATION
                        Path to GTF gene annotation file
  -d, --docker          Run all commands in the refgenie docker container.
  -o OUTFOLDER, --outfolder OUTFOLDER
                        Path to output genomes folder, using the $GENOMES
                        environment variableif set. Currently set to:
                        '/ext/yeti/genomes/'
```

