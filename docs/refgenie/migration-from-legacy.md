# Migration reference: 0.x to 1.0

This page lists the changes in refgenie 1.0 that affect code, pipelines, and plugins written for 0.x. To move an existing install, use the [upgrade guide](upgrade.md).

## Environment variables

Refgenie 1.0 does not read `$REFGENIE`. Use `$REFGENIE_HOME_PATH` or `$REFGENIE_DB_CONFIG_PATH` instead. While `$REFGENIE` is set and `$REFGENIE_DB_CONFIG_PATH` is not, refgenie prints a warning. See [Configure](configuration.md) for all variables.

## Seek keys

Every 0.x asset had an implicit `dir` seek key that returned its folder. In 1.0, each asset class declares typed seek keys (`file`, `directory`, `prefix`, `string`, `json`). Most asset classes have one main seek key, named like the class, which refgenie uses by default. So for most uses of `.dir`, drop the suffix:

| 0.x | 1.0 |
|---|---|
| `hg38/bowtie2_index.dir` (caller appends `/{genome}`) | `hg38/bowtie2_index` (returns the index prefix) |
| `hg38/bwa_index.dir` (caller appends `/{genome}.fa`) | `hg38/bwa_index` (returns the index prefix) |
| `hg38/fasta.dir` | `hg38/fasta` (returns the FASTA file path) |

In a looper pipeline interface (for example, PEPATAC):

```diff
- --genome-index { refgenie[sample.genome].bowtie2_index.dir }
+ --genome-index { refgenie[sample.genome].bowtie2_index }
```

If a tool needs a directory (for example, STAR's `--genomeDir`), its asset class declares a `directory` seek key.

## Looper integration

The pre-submit hook that fills the `refgenie` namespace moved:

```diff
 pre_submit:
   python_functions:
-  - refgenconf.looper_refgenie_populate
+  - refgenie.integrations.looper.populate
```

It reads the database from `var_templates.refgenie_db_config` in the pipeline interface when that is set, and from `$REFGENIE_DB_CONFIG_PATH` otherwise.

## Asset paths

Asset files are stored by digest, under `data/<genome_digest>/<asset_group>/<asset_digest>/`. Refgenie also builds a human-readable tree of symlinks under `alias/<alias>/<asset_group>/<asset_name>/`, with file names rewritten to the alias. `seek` returns the alias path by default; `--abs` (or `seek(..., abs_path=True)` in Python) returns the digest path. Both point to the same files.

```console
$ refgenie seek hg38/bwa_index:0.7.19
<genome_folder>/alias/hg38/bwa_index/0.7.19/hg38

$ refgenie seek --abs hg38/bwa_index:0.7.19
<genome_folder>/data/<genome_digest>/bwa_index/<asset_digest>/<genome_digest>
```

## Aligner-index recipes

Aligner-index asset classes such as `bowtie2_index` and `bwa_index` declare their main seek key as `type: prefix`, which resolves to `<asset_folder>/<genome>`, the value `bowtie2 -x` and `bwa mem` expect. If you maintain a recipe for an aligner index:

- Declare the seek key as `type: prefix` with `value: "{genome}"`.
- In the recipe command, use the aligner's prefix option (for example, `bwa index -p <prefix>`) so the output files match.

## Python API

- The `refgenconf.RefGenConf` object is replaced by `refgenie.Refgenie`. Its constructor takes `database_config_path: str | Path | None`. See [the Refgenie Python object](refgenie_object.md).
- `rgc.asset.seek(...)` returns a `str`.

## Build logs

Build logs and other build files (the pypiper log, commands file, profile, `stats.yaml`, and the completion flag) live in their own tree, `<genome_folder>/builds/<genome_alias>/<asset_group>/<asset_name>/`, outside the asset folder. Snakemake workflows should get the completion-flag path from `rgc.build.target_template()`. The flag file name is `{genome}_{group}__{asset}.flag`.

## Plugins

Plugins still register functions in `refgenie.hooks.<hook>` entry-point groups. The hooks, the function signature, and where settings live are different.

| 0.x hook | 1.0 |
| --- | --- |
| `post_update` | `post_update`: once, after any operation that changed local assets, aliases, or genomes |
| `pre_pull` / `post_pull` | same, once per asset pulled |
| `pre_update` | removed |
| `pre_tag` / `post_tag` | removed |
| `pre_list` / `post_list` | removed |
| (none) | `pre_build` / `post_build` |

`refgenie plugins` shows an entry point in a removed group as `unknown hook (never fires)`.

A hook function is called as `func(rg, event)`, where `rg` is the `Refgenie` instance and `event` is a `refgenie.plugins.HookEvent`. A function that takes one argument raises `TypeError`, which refgenie logs with the plugin's name before moving on.

0.x:

```python
def update_config(rgc):
    for genome, assets in rgc.list_assets_by_genome().items():
        for asset in assets:
            print(genome, asset, rgc.seek(genome, asset))
```

1.0:

```python
from refgenie.models import AssetRegistryPathComponents


def update_config(rg, event):
    paths = rg.paths()
    for genome in paths:  # genome aliases
        for group in paths[genome]:  # groups with a default asset
            path = rg.asset.seek_components(
                AssetRegistryPathComponents(genome=genome, asset_group=group)
            )
            print(genome, group, path)
```

Plugin settings live in the database, not in config-file keys (such as refgenie_nfcore's `nextflow_config:`). The user runs `refgenie plugins set nfcore config_path=/abs/path/nf.config`, and the plugin reads `rg.plugins.settings("nfcore")` in its hook. See [Write a plugin](plugins.md).
