# Refgenie configuration file upgrades

Refgenie is under active development and new features are added regularly. This sometimes necessitates changes in the refgenie configuration file format or asset directory structure.

Starting with the refgenie transition 0.9.3 -> 0.10.0 (configuration file versions: 0.3 -> 0.4) we introduced the `refgenie upgrade` functionality, which will take care of all the required reformatting. Running `refgenie upgrade` will automatically detect the current configuration file version and will: 1. reformat the configuration file to the new version; and 2) make any necessary changes to the asset directory structure.

Below we describe the changes introduced in each configuration file version and how to upgrade:

## Configuration file v0.4 (introduced: refgenie v0.10.0)

### How to upgrade

To reformat the config run from the command line:

```
refgenie upgrade --target-version 0.4 -c /path/to/old/cfg.yml
```

Or from within Python:

```python
from refgenconf import upgrade_config
upgrade_config(target_version="0.4", filepath="/path/to/old/cfg.yml")
```

### Config format changes

- use sequence-derived unique genome identifiers instead of genome names everywhere
- add `aliases` key under each genome section to store the aliases that can be used to refer to the genomes easily

### File tree structure changes

- use sequence-derived unique genome identifiers instead of genome names in every file name and directory name
- move all the contents from the refgenie directory to a new `data` directory
- add an `alias` directory with contents corresponding to the aliases defined in the configuration file. The contents of the child directories are symbolic links to the asset files in the `data` directory
