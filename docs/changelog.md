# Changelog

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html) and [Keep a Changelog](https://keepachangelog.com/en/1.0.0/) format. 

## [0.6.0] - 2019-08-05

### Added
- `list` and `listr` subcommand results can be restricted to a specific genome with `-g/--genome` options
- `remove` subcommand will remove an asset from disk and config
- Added recipes for new assets: `refgene_anno`, `ensembl_gtf` and `feat_annotation`
- `build` now populates the `asset_description` field in the config with corresponding value from the recipe


## [0.5.0] - 2019-07-11
### Changed
- `refgenie build` uses dict-like recipies for build instructions
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
