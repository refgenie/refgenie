# Changelog

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html) and [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) format.

## [0.13.0] - unreleased

### Added
- `--pipeline-kwargs` arguemnt to `refgenie build` command, which allows you to pass keyword arguments to the pypiper pipeline.

## [0.12.1] - 2021-11-04

### Fixed
- A bug with setuptools no longer allowing use_2to3.

## [0.12.0] - 2021-06-28

### Added
- _MapReduce_ framework support in `refgenie build` that supports automating asset builds at scale (`--map` and `--reduce` options)
- an option to automatically pull parent assets when building a derived asset (`--pull-parents` option in `refgenie build`)

### Fixed
- misleading exit codes in `refgenie build`; [#258](https://github.com/refgenie/refgenie/issues/258)

## [0.11.0] - 2021-04-27

### Added
- remote commands, which work without genome config file:
  - `refgenie seekr`
  - `refgenie populater`
  - `refgenie listr`
- `refgenie populate` command for refgenie registry paths populating with local paths
- `--flag-meanings` argument to `refgenie compare` command

## [0.10.0] - 2021-03-11

**Version `0.10.0` requires an upgrade to the configuration file and asset structure. Please refer to the [upgrade tutorial](config_upgrade_03_to_04.md) for instructions on how to migrate your config upon upgrade to `0.10.0`.**

### Changed

- instead of using human-readable names as genome identifiers refgenie uses sequence-derived digests in the config
- asset data moved to `data` directory
- asset files are now named after genome digests
- refgenieserver API v3 is now used for remote assets retrieval
- improved visual interface in `list`, `listr` and `pull` subcommands

### Added

- `data` and `alias` directories in genome directory that are used to store asset and aliases data, respectively
- `refgenie alias` command for genome aliases management
- `refgenie upgrade` command for config format upgrades
- `refgenie compare` command for genome compatibility determination

## [0.9.3] - 2020-07-29

### Changed
- short option string for `--no-overwrite` from `-n` to `-o`

### Added
- option to handle large asset archives pulling from the CLI (`-l`/`--no-large` flag)
- option to set the maximum archive size to `pull` with no confirmation required (`--size-cutoff` argument)
- `-s`/`--seek-keys` argument to `refgenie add` to specify seek keys for added assets

### Fixed
- `refgenie add` issues -- added assets are no longer imported to the `genome_folder`; [#180](https://github.com/refgenie/refgenie/issues/180)

## [0.9.2] - 2020-07-01

### Changed
- in `refgenie build` reduced the config file locking time to prevent problems in multi-build context
- dropped Python 2 support

### Added
- parametrized `kmer` in salmon recipes
- support for all genome configuration file parameter values initialization in `refgenie init`

## [0.9.1] - 2020-05-01

### Added
- added option (`-f`/`--force`) to confirm assets overwriting upfront in `refgenie add` add `refgenie pull`

### Changed
- fixed bug in hisat2_index that pointed to the parent folder. The seek key now points to the folder/{genome}, as expected by the tool
- fixed bug in bwa_index that pointed to the parent folder. The seek key now points to the folder/{genome}.fa, as expected by the tool

## [0.9.0] - 2020-03-17

### Changed
- fixed a bug in bowtie2_index recipe that pointed to the parent folder. The seek key now points to the folder/{genome}, as expected by bowtie2
- in `refgenie seek` file existence check is not performed by default

### Added
- possibility to execute library module as a script: `python -m refgenie ...`
- support for repeated recipe inputs on CLI (arguments: `--files`, `--assets` and `--params`)
- a possibility to perform file existence check (`-e`/`--check-exists`) in `refgenie seek`

## [0.8.2] - 2020-01-08

### Fixed
- `SyntaxError` in Python 2.7; [#155](https://github.com/databio/refgenie/issues/155)

## [0.8.1] - 2019-12-13

### Fixed
- `salmon_partial_sa_index` recipe

### Changed
- `refgenie remove` removes the asset relatives links
- `refgenie init` uses `initialize_config_file` method from `refgenconf`
- default input assets for `salmon_sa_index` and `salmon_partial_sa_index` recipes to transcriptomes within the namespace

### Added
- `threads` parameter to the following recipes: `dbnsfp`, `salmon_index`, `star_index`

### Removed
- documentation regarding `-r`/`--recipe` option in `refgenie build`. It will be removed in the future

## [0.8.0] - 2019-12-06

### Changed
- `refgenie build` command arguments naming scheme: `--{input_name} <path>` to `--files {input_name}=<path>`
- `-r`/`--requirements` in `refgenie build` command to `-q`/`--requirements`
- recipe format: requirements (both assets and inputs) are lists of dicts rather that lists of strings
- `refgenie list` displays current server subscriptions

### Added
- `refgenie id` command for asset digest retrieval
- cross-namespace asset relationships support
- `--assets` argument in `refgenie build` command to provide parent assets, if required
- `-r`/`--recipe` in `refgenie build` command argument to provide the recipe for the build
- `subscribe` and `unsubscribe` subcommands to enable server list manipulation in the config file (`genome_servers` entry in the refgenie configuration file)
- new recipes:
    - `salmon_sa_index`
    - `salmon_partial_sa_index`
    - `suffixerator_index`
    - `tallymer_index`

### Removed
- `-t`/`--tag` in `refgenie build`. Use more flexible `--assets` instead.

## [0.7.2] - 2019-11-06

### Added
- `dbsnp` recipe
- distribute the license file with the package


## [0.7.1] - 2019-10-29

### Changed
- `--genome-server` can now be called multiple times to add additional refgenieservers
- `listr` will check each available refgenieserver and display assets
- `pull` will check each available refgenieserver and take the first matching asset found

### Added
- possibility to list **multiple** selected genomes in `refgenie list/listr -g`

## [0.7.0] - 2019-10-21

### Added
- `import_igenome` command line tool for iGenomes integration with Refgenie
- `--genome/tag-description` arguments to the `refgenie build` command
- `-r`/`--requirements` argument to the `refgenie build` command recipe requirements to display required inputs and required assets for a particular recipe
- config manipulation support in multi-user contexts, it's racefree, uses file locks
- `dbNSFP` asset recipe
- assets tagging
- `refgenie tag` command that assigns a tag to an assets (re-tags it)
- `refgenie getseq` command that retrieves sequence ranges from a genome
- `seek_keys`, which provide control over files within an asset
- `asset_digests`, which are calculated after asset building and used to assure asset provenance
- asset relationships recording (`asset_children`, `asset_parents` fields)

### Changed
- assets can be referred to by registry paths: `genome/asset.seek_key:tag`
- config v0.3 is required
- `refgenie pull` uses `refgenieserver` API v2

## [0.6.0] - 2019-08-05

### Added
- `list` and `listr` subcommand results can be restricted to a specific genome with `-g/--genome` options
- `remove` subcommand will remove an asset from disk and config
- Added recipes for new assets: `ensembl_gtf` and `feat_annotation`
- `build` now populates the `asset_description` field in the config with corresponding value from the recipe

### Changed
- changed some asset locations; `tss_annotation` is now named `refgene_tss` or `ensembl_tss` and is built by the `refgene_anno` or `ensembl_gtf` assets. Renamed `gene_anno` to `refgene_anno`.


## [0.5.0] - 2019-07-11
### Changed
- `refgenie build` uses dict-like recipes for build instructions
- Major genome configuration file format changes
    - Added `config_version` entry
    - Added `assets` section in `genomes` section

- recipes can now include container images

### Added
- `genomes` can have attributes, like description
- Added recipes for new assets `salmon`, `bwa`, `star`, `gene_anno`, and `tss_annotation`.

## [0.4.4] - 2019-07-01
### Added
- `add` subcommand

## [0.4.3] - 2019-06-21
### Changed
- Build process now builds individual assets

## [0.4.2] - 2019-06-18
### Added
- `seek` subcommand

### Changed
- Require config file arg for `init`.

### Fixed
- Pick up env var for `init` config.

## [0.4.1] -- 2019-06-14
### Fixed
- Use newer `yacman` and regain `init` functionality.

## [0.4.0] -- 2019-06-14
### Added
- Added new commands `init`, `pull`, `list` and `listr`
- Added connectivity option with remote data sources

## [0.3.2] -- 2019-05-14
### Fixed
- Fixed a bug with packaging

## [0.3.0] -- 2019-05-10
### Added
- Implemented installable CLI
- Packaged for release on PyPI
- Wrote comprehensive docs
### Fixed
- Fixed naming of `.fq.gz` files

## [0.2.0] -- 2017-03-08
### Added
- Transition release as a functional script
## [0.1.0] -- 2016-11-11
- Project started
