# Use refgenie in your pipeline

The code snippets below can be used in your pipeline to **assert the existence of the refgenie-managed files** in 3 different languages: Bash, Python and R.

Refgenie checks if the asset is available locally and tries pull it from the server if it's not.

_The only_ step that needs to precede the execution of these functions is refgenie genome configuration file initialization:

```console
export REFGNEIE=refgenie_config.yaml
refgenie init -c $REFGENIE
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
    if [ -z "$REFGENIE" ]
    then
        echo -e "${RED}refgenie env var not defined."
        echo -e "Run 'export REFGENIE=<path to refgenie config>' to set the env var.${NC}"
        exit 1
    else
        echo -e "${GREEN}refgenie env var defined: $REFGENIE${NC}"
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

* Python package `refgenconf`

```python
from refgenconf import RefGenConf

def assert_refgenie_asset_exists(
    genome, asset, tag=None, seek_key=None, refgenie_config=None
):
    # instantiate RefGenConf object
    rgc = RefGenConf(filepath=refgenie_config)

    # get tag of interest, provided vs. default
    tag = tag if tag is not None else rgc.get_default_tag(genome=genome, asset=asset)

    # list assets available locally
    list_result = rgc.list()

    # check whether the asset of interest is missing
    if genome not in list_result.keys() or asset not in list_result[genome]:
        # pull asset if missing
        print(f"{genome}/{asset}:{tag} not found, pulling...")
        try:
            rgc.pull(genome=genome, asset=asset, tag=tag)
        except Exception as e:
            print(f"Pull failed")
            raise

    # get the local path to the asset of interest
    rgc.seek(genome=genome, asset=asset, tag=tag, seek_key=seek_key)


# Run like this: assert_refgenie_asset_exists(
#     genome="hg38",
#     asset="fasta",
# )
```

## R

**Requirements:**

* Python package `refgenconf`
* R package `reticulate`


```
library('reticulate')

assertRefgenieAssetExists <-
  function(genome,
           asset,
           tag = NULL,
           seek_key = NULL,
           refgenieConfig = NULL) {

    # import Python module
    refgenconf = reticulate::import("refgenconf", convert = FALSE)

    # determine refgenie config path, provided vs. read from env
    refgenieConfig = ifelse(is.null(refgenieConfig),
                            Sys.getenv("REFGENIE"),
                            refgenieConfig)

    # instantiate Python RefGenConf object
    rgc = refgenconf$RefGenConf(filepath = refgenieConfig)

    # get tag of interest, provided vs. default
    tag = ifelse(is.null(tag),
                 py_to_r(rgc$get_default_tag(genome = genome, asset = asset)),
                 tag)

    # string together the final asset registry path, for logging
    assetRegistryPath = paste0(genome, "/" , asset, ":",  tag)

    # list assets available locally
    listResult = py_to_r(rgc$list())

    # check whether the asset of interest is missing
    if (is.null(listResult[[genome]]) |
        !any(listResult[[genome]] == asset)) {
      # pull asset if missing
      message(paste0(assetRegistryPath, " not found, pulling..."))
      pullResult = py_to_r(rgc$pull(
        genome = genome,
        asset = asset,
        tag = tag,
        force = TRUE,
        force_large = TRUE
      ))
    }

    # get the local path to the asset of interest
    seekResult = rgc$seek(
      genome_name = genome,
      asset_name = asset,
      tag_name = tag,
      seek_key = seek_key
    )
  }


# Run like this: assertRefgenieAssetExists(
#     genome="hg38",
#     asset="fasta",
# )
```
