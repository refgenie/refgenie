# Manage local assets

The `init` command initializes an empty genome config file for you:

```
refgenie init -c genome_config.yaml
```


The `list` command *lists local assets*:

```
refgenie list
```

## The REFGENIE environment variable

The `refgenie` commands all require knowing where your genome configuration file is. You can pass it on the command line all the time (using the `-c` parameter), but this gets old if you use refgenie a lot. An alternative is to simply set up the $REFGENIE environment variable like so:

```
export REFGENIE=/path/to/genome_config.yaml
```

Add this to your `.bashrc` or `.profile` if you want it to persist. This way, you won't have to specify `-c` to every call. You can always specific `-c` if you want to override the value in the $REFGENIE variable.