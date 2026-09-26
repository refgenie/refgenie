# Data Channels Configuration

## Overview

Refgenieserver supports data channels - static files hosted by third-party providers that contain asset classes and recipes. Data channels are aggregated into a single unified interface, making it transparent to end users.

## Configuration File

Create a `data_channels.yaml` file to configure multiple data channels:

```yaml
# data_channels.yaml
main:
  name: main
  protocol: https
  index_url: https://refgenie.github.io/recipes/index.yaml

custom_channel:
  name: custom_channel
  protocol: https
  index_url: https://example.com/my-recipes/index.yaml

dev:
  name: dev
  protocol: http
  index_url: http://localhost:8080/index.yaml
```

## Configuration Parameters

- **name**: Channel identifier (cannot contain `__` characters)
- **protocol**: `http` or `https`
- **index_url**: URL to the channel's `index.yaml` file

## Environment Configuration

Set the configuration file path via environment variable:

```bash
export DATA_CHANNELS_CONFIG_PATH="/path/to/your/data_channels.yaml"
```

Default location: `refgenieserver/data_channels.yaml`

## Data Channel Format

Each data channel must expose an `index.yaml` file with this structure:

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

## API Endpoints

- `GET /data_channels/` - List all configured channels
- `GET /data_channels/index.yaml` - Compiled index from all channels
- `GET /data_channels/{channel_name}__{filename}` - Access files from specific channels

## How It Works

1. Refgenieserver loads all configured data channels at startup
2. Creates an aggregated index combining all channel files
3. Files are accessed using the format: `{channel_name}__{filename}`
4. Requests are transparently redirected to the appropriate data channel URL