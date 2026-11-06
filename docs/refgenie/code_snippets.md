# Use refgenie in your pipeline

These snippets make sure a refgenie asset is present before a pipeline uses it, in Bash, Python, and R. Each one looks for the asset locally, pulls it from your server if it is missing, and returns its path.

Run `refgenie init` once first. By default refgenie uses a local SQLite database; to use another database, set `REFGENIE_DB_CONFIG_PATH` to its configuration file.

## Bash

Needs the `refgenie` Python package.

```bash
#!/bin/bash

assert_refgenie_asset_exists(){
    if [ -z "$1" ]; then
        echo "Asset registry path not provided!" >&2
        exit 1
    fi

    if file_path=$(refgenie seek "$1"); then
        echo "Found ($1) asset: $file_path" >&2
    else
        echo "Asset ($1) not found, pulling..." >&2
        refgenie pull "$1"
        if ! file_path=$(refgenie seek "$1"); then
            echo "Asset ($1) pull failed" >&2
            exit 1
        fi
    fi
    echo "$file_path"
}

# Run like this: assert_refgenie_asset_exists hg38/fasta
```

## Python

Needs the `refgenie` Python package.

```python
from refgenie import GenomeAlias, Refgenie
from refgenie.exceptions import MissingAliasError


def assert_refgenie_asset_exists(
    genome, asset_group, asset=None, seek_key=None, refgenie_config=None
):
    # uses the default database when refgenie_config is None
    rg = Refgenie(database_config_path=refgenie_config)
    alias = GenomeAlias(genome)

    try:
        digest = rg.alias.resolve(alias)
        have_group = rg.asset.group.exists(asset_group, genome_digest=digest)
    except MissingAliasError:
        have_group = False

    if not have_group:
        print(f"{genome}/{asset_group} not found, pulling...")
        rg.transfer.pull(asset_group_name=asset_group, genome=alias, asset_name=asset)
        digest = rg.alias.resolve(alias)

    # the group's default asset and seek key unless given
    return rg.asset.seek(digest, asset_group, asset_name=asset, seek_key_name=seek_key)


# Run like this: assert_refgenie_asset_exists(genome="hg38", asset_group="fasta")
```

## R

Needs the `refgenie` Python package and the R package [reticulate](https://rstudio.github.io/reticulate/).

```r
library(reticulate)

assertRefgenieAssetExists <- function(genome, assetGroup, asset = NULL, seekKey = NULL) {
  refgenie <- reticulate::import("refgenie")

  rg <- refgenie$Refgenie()
  alias <- refgenie$GenomeAlias(genome)

  haveGroup <- tryCatch({
    digest <- rg$alias$resolve(alias)
    rg$asset$group$exists(assetGroup, genome_digest = digest)
  }, error = function(e) FALSE)

  if (!haveGroup) {
    message(paste0(genome, "/", assetGroup, " not found, pulling..."))
    rg$transfer$pull(asset_group_name = assetGroup, genome = alias, asset_name = asset)
  }

  digest <- rg$alias$resolve(alias)
  rg$asset$seek(digest, assetGroup, asset_name = asset, seek_key_name = seekKey)
}

# Run like this: assertRefgenieAssetExists(genome = "hg38", assetGroup = "fasta")
```
