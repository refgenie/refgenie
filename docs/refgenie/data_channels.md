# Refgenie Data Channels

## Introduction

Refgenie data channels are a powerful mechanism for extending and sharing asset classes and recipes across the refgenie ecosystem. Data channels enable tool developers, organizations, and the community to publish, distribute, and synchronize new types of reference genome assets and build instructions without modifying the core refgenie codebase.

!!! success "Learning objectives"
    - What are refgenie data channels?
    - Why are data channels useful?
    - How do I connect my local refgenie to a remote data channel?
    

!!! warning "Warning"
    This guide only shows you how to *use existing channels*. If you want to set up your own channel, please see [Set up your own data channel](set_up_data_channel.md)


## What is a data channel?

A data channel is an external source (typically a remote repository or URL) that provides definitions for asset classes and recipes. These definitions describe how to build, manage, and use new types of reference genome assets. Data channels are registered with your local refgenie instance and can be synchronized to keep your asset classes and recipes up to date.

- **Asset class**: Defines the structure and seek keys for a type of asset (e.g., `fasta`, `bwa_index`).
- **Recipe**: Describes how to build an asset class from input files and parameters.

A data channel is basically a set of `.yaml` files representing asset classes and recipes, listed in an `index.yaml` file, like this:

```yaml
asset_class:
  dir: asset_classes
  files:
    - fasta.yaml
    - bowtie2_index.yaml
recipe:
  dir: recipes
  files:
    - fasta.yaml
    - bowtie2_index.yaml
```

## Why use data channels?

- **Community-driven extension**: Anyone can publish new asset classes and recipes for others to use.
- **Separation of concerns**: Core refgenie remains stable while new types of assets and related recipes are distributed via channels.
- **Easy updates**: Syncing a channel brings in the latest definitions without reinstalling refgenie.
- **Reproducibility**: Standardized recipes and asset classes can be shared and reused across projects and organizations.

## Example data channel

Refgenie authors maintain a demo data channel at [https://refgenie.github.io/recipes](https://refgenie.github.io/recipes). This channel is served via GitHub Pages and provides a curated collection of asset class and recipe definitions, an index file, and a simple landing page. It is designed as a reference implementation and a starting point for anyone wishing to distribute their own asset classes and recipes.

**Key features of the default data channel:**

- **Simple structure:** The channel consists of YAML files for asset classes and recipes, organized and indexed for easy consumption by refgenie.
- **Ready-to-use:** By default, this channel provides definitions for common assets such as `fasta`, `bwa_index`, and more, making it easy to get started.
- **Template for custom channels:** Tool developers and users are encouraged to fork or mimic this structure to create their own data channels for organizational or specialized needs.

**How to use:**

- You can register the demo channel directly with:

  ```bash
  refgenie data-channel add demo https://refgenie.github.io/recipes/index.yaml
  refgenie data-channel sync demo --exists-ok
  ```

- After syncing, all asset classes and recipes from the channel will be available for use in your local refgenie instance.

For more information or to contribute, visit the [GitHub repository](https://github.com/refgenie/recipes).

## Registering and Syncing Data Channels

You can add a data channel to your refgenie instance using the CLI:

```bash
refgenie data-channel add my-channel \
  https://refgenie.github.io/recipes/index.yaml
```

List registered data channels:

```bash
refgenie data-channel list
```

To synchronize (download and register) all asset classes and recipes from a channel:

```bash
refgenie data-channel sync my-channel --exists-ok
```

This will fetch the latest definitions and make them available for building and managing assets.

## Using Data Channels Programmatically

You can also interact with data channels from Python:

```python
from refgenie import Refgenie
refgenie = Refgenie()

# Add a data channel
refgenie.sources.add_channel(
    name="my-channel",
    type="http",
    index_address="https://refgenie.github.io/recipes/index.yaml",
    description="Community recipes"
)

# List asset classes and recipes from the channel
for asset_class in refgenie.sources.iter_asset_classes("my-channel"):
    refgenie.asset_class.add(asset_class)
for recipe in refgenie.sources.iter_recipes("my-channel"):
    refgenie.recipe.add(recipe)
```

---

## Data Channels in Refgenieserver

Refgenieserver is a companion server software for refgenie that can host and aggregate data channels for your organization or community. When you point your refgenie client to a refgenieserver instance, you gain access to all asset classes and recipes provided by the server's configured data channels—without having to manage multiple URLs or repositories yourself.

**Key points for refgenie users:**

- You can use a single refgenieserver endpoint to access asset classes and recipes from many data channels, simplifying configuration and collaboration.
- The server automatically keeps its index up to date, so you always see the latest available definitions.
- From the user's perspective, working with data channels via refgenieserver is seamless—just subscribe to the server and use the available assets and recipes as usual.

This makes refgenieserver ideal for organizations, consortia, or communities that want to provide a unified, up-to-date set of reference genome resources to all their users.

