# How to use `refgenconf ` to manage Refgenie assets in a pipeline


Below we present an example use of `refgenconf` package. It is installed automatically with `refgenie` (or separately installable with `pip install refgenconf`). All the asset fetching functionality is impelmented in `refgenconf` package, so pipelines that just use Python API do not need to depend on `refgenie`.

## Goal
The goal of the code below is to **get a path to the refgenie-managed fasta file for a user-specified genome**. 

Genome FASTA is a part of `fasta` asset, accessible as a `fasta` seek key. To retrieve the path this file on the command line one would say: `refgenie seek <genome>/fasta`. For example:
```
refgenie seek hg38/fasta
```

## Steps

First, let's set the `$REFGENIE` environmet variable. It should be set by a pipeline user or the config file path should be provided explictly, e.g. as an input to the pipeline (here shown as `user_provided_cfg_path = None` -- not provided) 


```python
import os
os.environ["REFGENIE"] = "./refgenie.yaml"
user_provided_cfg_path = None
user_provided_genome = "rCRSd"
```

Next, let's import components of `refgenconf` that we'll use


```python
from refgenconf import RefGenConf, select_genome_config, RefgenconfError, CFG_ENV_VARS, CFG_FOLDER_KEY
from yacman import UndefinedAliasError
```

Now, we can use the `select_genome_config` function to determine the preferred path to the config file. If `user_provided_cfg_path` is `None` (not specified) the `$REFGENIE` environment variable is used. 


```python
refgenie_cfg_path = select_genome_config(filename=user_provided_cfg_path, check_exist=False)
```

The function returns `None` if none of the above point to a valid path. That's why we raise an aproppriate error below. Obviously, the name of `--rfg-config` argument depends on pipeline design. 


```python
if not refgenie_cfg_path:
    raise OSError(f"Could not determine path to a refgenie genome configuration file."
                  f"Use --rfg-config argument or set '{CFG_ENV_VARS}' environment variable to provide it")
```

Otherwise it returns a determined path (`str`). So, we check if it exists and read the object if it does. If it does not, we can initialize the config file


```python
if isinstance(refgenie_cfg_path, str) and os.path.exists(refgenie_cfg_path):
    print(f"Reading refgenie genome configuration file from file: {refgenie_cfg_path}")
    rgc = RefGenConf(filepath=refgenie_cfg_path)
else:
    print(f"File '{refgenie_cfg_path}' does not exist. Initializing refgenie genome configuration file.")
    rgc = RefGenConf(entries={CFG_FOLDER_KEY: os.path.dirname(refgenie_cfg_path)})
    rgc.initialize_config_file(filepath=refgenie_cfg_path)
    rgc.subscribe(urls="http://rg.databio.org:82", reset=True) # subscribe to the desired server, if needed
```

    File '/Users/mstolarczyk/code/refgenie/docs_jupyter/refgenie.yaml' does not exist. Initializing refgenie genome configuration file.


Finally, we try to retrieve the path to out asset of interest and pull from `refgenieserver` if the retrieval fails.


```python
try:
    fasta = rgc.seek(genome_name=user_provided_genome, asset_name="fasta", tag_name="default",
                                seek_key="fasta")
except (RefgenconfError, UndefinedAliasError):
    print("Could not determine path to chrom.sizes asset, pulling")
    rgc.pull(genome=user_provided_genome, asset="fasta", tag="default")
    fasta = rgc.seek(genome_name=user_provided_genome, asset_name="fasta", tag_name="default",
                                seek_key="fasta")
print(f"Determined path to fasta asset: {fasta}")
```

    Could not determine path to chrom.sizes asset, pulling



    Output()


    Determined path to fasta asset: /Users/mstolarczyk/code/refgenie/docs_jupyter/alias/rCRSd/fasta/default/rCRSd.fa



```python

```
