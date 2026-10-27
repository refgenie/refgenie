# Use refgenie in your pipeline

The code snippets below can be used in your pipeline to **assert the existence of the refgenie-managed files** in 3 different languages: Bash, Python and R.

Refgenie checks if the asset is available locally and tries pull it from the server if it's not.

_The only_ step that needs to precede the execution of these functions is initializing refgenie's database configuration:

```console
refgenie init
```

By default this creates a local SQLite database. To use a specific configuration file, set the `REFGENIE_DB_CONFIG_PATH` environment variable to its path:

```console
export REFGENIE_DB_CONFIG_PATH=/path/to/refgenie_db_config.yaml
```

## Bash

**Requirements:**

* Python package `refgenie`


```bash
#!/bin/bash

assert_refgenie_asset_exists(){

    RED='\033[0;31m'
    GREEN='\033[0;32m'
    YELLOW='\033[0;33m'
    NC='\033[0m'

    if [ -z "$1" ]; then
        echo -e "\n${RED}Asset registry path not provided!${NC}\n"
        exit 1
    fi

    # check if refgenie env var is defined
    if [ -z "$REFGENIE_DB_CONFIG_PATH" ]
    then
        echo -e "${RED}refgenie env var not defined."
        echo -e "Run 'export REFGENIE_DB_CONFIG_PATH=<path to refgenie db config>' to set the env var.${NC}"
        exit 1
    else
        echo -e "${GREEN}refgenie env var defined: $REFGENIE_DB_CONFIG_PATH${NC}"
    fi

    # check if asset is available locally
    if file_path=`refgenie seek $1`; then
        echo -e "${GREEN}Found ($1) asset: $file_path${NC}"
    else
        # pull if not available locally
        echo -e "${YELLOW}Asset ($1) not found, pulling...${NC}"
        refgenie pull $1
        if file_path=`refgenie seek $1`; then
            echo -e "${GREEN}Asset ($1) pulled successfully: $file_path${NC}"
        else
            echo -e "${RED}Asset ($1) pull failed${NC}"
            exit 1
        fi
    fi
}


# Run like this: assert_refgenie_asset_exists hg38/fasta
```


## Python

**Requirements:**

* Python package `refgenie`

```python
from refgenie import GenomeAlias, Refgenie
from refgenie.exceptions import MissingAliasError

def assert_refgenie_asset_exists(
    genome, asset_group, asset=None, seek_key=None, refgenie_config=None
):
    # instantiate Refgenie object (defaults to the local database if config is None)
    rg = Refgenie(database_config_path=refgenie_config)
    alias = GenomeAlias(genome)

    # check whether the genome or the asset group is missing locally
    try:
        digest = rg.alias.resolve(alias)
        have_group = rg.asset.group.exists(asset_group, genome_digest=digest)
    except MissingAliasError:
        have_group = False

    if not have_group:
        # pull asset if missing
        print(f"{genome}/{asset_group} not found, pulling...")
        try:
            rg.pull(asset_group_name=asset_group, genome=alias, asset_name=asset)
        except Exception as e:
            print("Pull failed")
            raise
        digest = rg.alias.resolve(alias)

    # get the local path to the asset of interest (the group's default asset if none given)
    return rg.asset.seek(
        digest,
        asset_group,
        asset_name=asset,
        seek_key_name=seek_key,
    )


# Run like this: assert_refgenie_asset_exists(
#     genome="hg38",
#     asset_group="fasta",
# )
```

## R

**Requirements:**

* Python package `refgenie`
* R package `reticulate`


```
library('reticulate')

assertRefgenieAssetExists <-
  function(genome,
           assetGroup,
           asset = NULL,
           seek_key = NULL,
           refgenieConfig = NULL) {

    # import Python module
    refgenie = reticulate::import("refgenie", convert = FALSE)

    # determine refgenie db config path, provided vs. read from env
    refgenieConfig = ifelse(is.null(refgenieConfig),
                            Sys.getenv("REFGENIE_DB_CONFIG_PATH"),
                            refgenieConfig)

    # instantiate Python Refgenie object
    rgc = refgenie$Refgenie(database_config_path = refgenieConfig)

    # get the asset (tag) of interest, provided vs. default
    asset = ifelse(is.null(asset),
                   py_to_r(rgc$asset$get_default(assetGroup, genome_name = genome)),
                   asset)

    # string together the final asset registry path, for logging
    assetRegistryPath = paste0(genome, "/" , assetGroup, ":",  asset)

    # check whether the asset group is missing locally
    if (!py_to_r(rgc$asset$group_exists(assetGroup, genome_name = genome))) {
      # pull asset if missing
      message(paste0(assetRegistryPath, " not found, pulling..."))
      pullResult = py_to_r(rgc$pull(
        asset_group_name = assetGroup,
        alias_name = genome,
        asset_name = asset
      ))
    }

    # get the local path to the asset of interest
    seekResult = rgc$asset$seek(
      genome_name = genome,
      asset_group_name = assetGroup,
      asset_name = asset,
      seek_key_name = seek_key
    )
  }


# Run like this: assertRefgenieAssetExists(
#     genome="hg38",
#     assetGroup="fasta",
# )
```
