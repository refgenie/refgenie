# Data Channels in Refgenieserver

## Overview

Data channels provide a way for a third-party user to provide recipes and asset classes to the refgenie ecosystem.
Refgenie users or servers can then subscribe to different channels, which makes available the recipes from those channels.




## Configuring Data Channels

Data channels in refgenieserver are configured using a YAML file that lists all available channels and their properties. When deploying a refgenieserver instance, you can specify the path to this configuration file using the `DATA_CHANNELS_CONFIG_PATH` environment variable. If this variable is not set, the server will look for a default file (typically named `data_channels.yaml`) in the project directory.

At server startup, the configuration is loaded and validated. If the configuration file is missing, empty, or invalid, the server will not start and will raise an error. This ensures that only valid data channel configurations are used.

The server also provides an endpoint to list all available data channels, making it easy for users and administrators to verify which channels are currently active.

## How Refgenieserver Handles Data Channels

- **Static file hosting:** Data channels are simply static files (YAML definitions) hosted by third-party providers. Each channel must expose an `index.yaml` file listing all available asset class and recipe files, following the standard format:

  ```yaml
  asset_class:
    dir: asset_classes
    files:
      - bowtie2_index_asset_class.yaml
      - bwa_index_asset_class.yaml
      - fasta_asset_class.yaml
  recipe:
    dir: recipes
    files:
      - bowtie2_index_asset_recipe.yaml
      - bwa_index_asset_recipe.yaml
      - fasta_asset_recipe.yaml
  ```

- **Channel aggregation:** Refgenieserver can be configured with multiple data channels (a dictionary mapping channel names to their `index.yaml` URLs). At server startup, it compiles an aggregated in-memory index from all configured channels, making all definitions available through a unified interface.

- **Transparent routing:** When a user requests a file (e.g., an asset class or recipe), refgenieserver transparently redirects the request to the appropriate data channel and file. To the end user, it appears as if they are interacting with a single, seamless data channel, even if the content is aggregated from multiple sources.

- **Dynamic updates:** The server can refresh and recompile the index as needed, ensuring users always have access to the latest definitions from all configured channels.

## Example Use Case

Suppose you have several organizations or teams each maintaining their own data channel. Refgenieserver can aggregate all these channels, so users only need to point their refgenie client to a single server endpoint to access all available asset classes and recipes.

## Benefits

- **Centralized access:** Aggregate and serve multiple data channels from a single endpoint.
- **Scalability:** Easily add or remove channels as your ecosystem grows.
- **Transparency:** Users do not need to know the details of each channel; routing is handled automatically.