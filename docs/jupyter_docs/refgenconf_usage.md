jupyter:True
#  How to `refgenconf ` to manage Refgenie assets in a pipeline

Below we present an example use of `refgenconf` package that is installed automatically with `refgenie` (or separately installable with `pip install refgenconf`). All the asset fetching functionality is impelmented in `refgenconf` package, so pipelines that just use Python API do not need to depend on `refgenie`.

## Goal
The goal of the code below is to **get a path to the refgenie-managed fasta file for a user-specified genome**. 

Genome FASTA is a part of `fasta` asset, accessible as a `fasta` seek key. To retrieve the path this file on the command line one would say: `refgenie seek <genome>/fasta`. For example:
```
refgenie seek hg38/fasta
```

## Steps

First, let's set the `$REFGENIE` environmet variable. It should be set by a pipeline user or the config file path shoulg be provided explictly, e.g. as an input to the pipeline (here shown as `user_provided_cfg_path`, not provided) 


```python
import os
os.environ["REFGENIE"] = "./genomes.yaml"
user_provided_cfg_path = None
user_provided_genome = "rCRSd"
```

Next, let's import components of `refgenconf` that we'll use


```python
from refgenconf import RefGenConf, select_genome_config, RefgenconfError, CFG_ENV_VARS, CFG_FOLDER_KEY
```

Now, we can use the `select_genome_config` function to determine the preferred path to the config file. If `user_provided_cfg_path` is `None` (not specified) the `$REFGENIE` environment variable is used. 


```python
refgenie_cfg_path = select_genome_config(filename=user_provided_cfg_path, check_exist=False)
```

The function returns `None` none of the above point to a valid path. That's why we raise an aproppriate error below 


```python
if not refgenie_cfg_path:
    raise OSError("Could not determine path to a refgenie genome configuration file. "
                  "Use --rfg-config argument or set '{}' environment variable to provide it".
                  format(CFG_ENV_VARS))
```

Otherwise it returns a determined path (`str`). So, we check if it exists and read the object if it does. If it does not, we can initialize the config file


```python
if isinstance(refgenie_cfg_path, str) and os.path.exists(refgenie_cfg_path):
    print("Reading refgenie genome configuration file from file: {}".format(refgenie_cfg_path))
    rgc = RefGenConf(filepath=refgenie_cfg_path)
else:
    print("File '{}' does not exist. Initializing refgenie genome configuration file.".format(refgenie_cfg_path))
    rgc = RefGenConf(entries={CFG_FOLDER_KEY: os.path.dirname(refgenie_cfg_path)})
    rgc.initialize_config_file(filepath=refgenie_cfg_path)
```

```.output
File '/Users/mstolarczyk/Uczelnia/UVA/code/refgenie/docs_jupyter/genomes.yaml' does not exist. Initializing refgenie genome configuration file.

```

Finally, we try to retrieve the path to out asset of interest and pull from `refgenieserver` if the retrieval fails.


```python
try:
    fasta = rgc.get_asset(genome_name=user_provided_genome, asset_name="fasta", tag_name="default",
                                seek_key="fasta")
except RefgenconfError:
    print("Could not determine path to chrom.sizes asset, pulling")
    rgc.pull_asset(genome=user_provided_genome, asset="fasta", tag="default")
    fasta = rgc.get_asset(genome_name=user_provided_genome, asset_name="fasta", tag_name="default",
                                seek_key="fasta")
print("Determined path to fasta asset: {}".format(fasta))
```

```.output
                                                      
```

```.output
Could not determine path to chrom.sizes asset, pulling
Determined path to fasta asset: /Users/mstolarczyk/Uczelnia/UVA/code/refgenie/docs_jupyter/rCRSd/fasta/default/rCRSd.fa

```

```.output

```
