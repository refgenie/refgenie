# Changelog

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html) and [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) format. 

## [0.7.2] - 2019-11-XX

### Added
- `dbsnp` recipe

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
