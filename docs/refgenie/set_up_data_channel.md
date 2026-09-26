# Setting Up Data Channels for Refgenie

## Introduction

This guide walks you through creating your own data channel to distribute custom asset classes and recipes for the refgenie ecosystem.

!!! success "Learning objectives"
    - How do I set up my own data channel?
    - How can I test a local or remote data channel to make sure it's set up correctly?


## Overview

A data channel is a collection of asset classes and recipes hosted at a URL that can be consumed by refgenie clients. Setting up your own data channel allows you to:

- Share custom asset classes and recipes with your organization or the community
- Maintain control over your genomic asset definitions
- Create specialized workflows for your specific use cases
- Contribute to the refgenie ecosystem

## Prerequisites

- A GitHub account (for GitHub Pages hosting) or web server
- Basic knowledge of YAML syntax
- Python 3.x installed (for validation tools)

## Step 1: Repository Structure

The easiest way to do this is to clone the [Official Recipes Repository](https://github.com/refgenie/recipes), and then just delete the recipes and asset classes and add your own. Or, you can create your own repository with the following structure:

```
my-data-channel/
    asset_classes/         # Directory for asset class definitions
        *.yaml             # Asset class YAML files
    recipes/               # Directory for recipe definitions
        *.yaml             # Recipe YAML files
    index.yaml             # Index file listing all available files
    build_index.py         # Script to generate index.yaml
    data_channel_check.py  # Script to validate channel content
    index.html             # Optional landing page
```

## Step 2: Create Asset Classes

Asset classes define the structure and seek keys for a type of asset. Create YAML files in the `asset_classes/` directory:

**Example: `asset_classes/my_index_asset_class.yaml`**

```yaml
name: my_index
version: 0.0.1
description: Custom index format for my tool
seek_keys:
  index_file:
    value: "{genome}.idx"
    description: Main index file
    type: file
  metadata:
    value: "{genome}.meta"
    description: Index metadata
    type: file
parents: []
```

For more details, see [Asset Class specification](asset_class_specification.md).

## Step 3: Create Recipes

Recipes describe how to build an asset class from inputs. Create YAML files in the `recipes/` directory:

**Example: `recipes/my_index_asset_recipe.yaml`**

```yaml
name: my_index
version: 0.0.1
output_asset_class: my_index
description: Build custom index from FASTA file
input_files:
  fasta:
    description: Input FASTA file (gzipped or not)
input_params: null
input_assets: null
docker_image: docker.io/myorg/my-indexer:latest
command_templates:
  - my-indexer build {{values.files["fasta"]}} -o {{values.output_folder}}/{{values.genome_digest}}.idx
  - my-indexer meta {{values.output_folder}}/{{values.genome_digest}}.idx > {{values.output_folder}}/{{values.genome_digest}}.meta
custom_properties:
  version: "my-indexer --version | head -1"
default_asset: "my-indexer-{{values.custom_properties.version}}"
```

For more details, see [Recipe specification](recipe_specification.md).

## Step 4: Generate the Index File

Next we need to create the index file. The `build_index.py` script will do this, and it will run automatically via a GitHub Action. If you like, you may also run it manually to generate your index:

```bash
python build_index.py
```

The generated `index.yaml` will look like:

```yaml
asset_class:
  dir: asset_classes
  files:
    - my_index_asset_class.yaml
recipe:
  dir: recipes
  files:
    - my_index_asset_recipe.yaml
```

## Step 5: Validate Your Data Channel

Use the included validation script to ensure your channel is properly formatted:

```bash
python data_channel_check.py .
```

## Step 6: Host Your Data Channel

### Option A: GitHub Pages (Recommended)

1. Push your repository to GitHub
2. Enable GitHub Pages in repository settings (Settings > Pages)
3. Select source branch (usually `main` or `master`)
4. Your channel will be available at: `https://[username].github.io/[repository-name]/`

You can test it with `python data_channel_check.py https://[username].github.io/[repository-name]/`

The demo repository also has GitHub Actions already set up to automatically update `index.yaml` when files change.

### Option B: Static Web Server

Upload your files to any static web server ensuring:
- All YAML files are accessible via HTTP/HTTPS
- CORS headers allow cross-origin requests (if needed)
- The `index.yaml` file is at the root of your channel URL


## Data Channel API Specification

Instead of hosting static files, you can implement a server that provides the Data Channel API. This allows for dynamic generation of asset classes and recipes, database-backed storage, authentication, and other advanced features.

### API Overview

A data channel server must provide the following HTTP endpoints:

#### 1. Index Endpoint

**URL:** `/index.yaml` (or base URL)
**Method:** GET
**Content-Type:** `text/yaml` or `application/x-yaml`

**Response Format:**
```yaml
asset_class:
  dir: asset_classes
  files:
    - file1.yaml
    - file2.yaml
recipe:
  dir: recipes
  files:
    - recipe1.yaml
    - recipe2.yaml
```

The index specifies:
- `dir`: The relative path to the directory containing the files
- `files`: List of available YAML files

#### 2. Asset Class Endpoints

**URL Pattern:** `/{asset_class_dir}/{filename}`
**Method:** GET
**Content-Type:** `text/yaml` or `application/x-yaml`

**Response:** Asset class definition YAML conforming to the [Asset Class specification](asset_class_specification.md)

#### 3. Recipe Endpoints

**URL Pattern:** `/{recipe_dir}/{filename}`
**Method:** GET
**Content-Type:** `text/yaml` or `application/x-yaml`

**Response:** Recipe definition YAML conforming to the [Recipe specification](recipe_specification.md)

### API Requirements

To be compatible with refgenie clients, your server must return valid YAML. All responses must be parseable YAML.


## Step 7: Use your data channel

Once hosted, test your channel with refgenie:

```bash
# Add your channel
refgenie data-channel add my-channel https://[username].github.io/[repository-name]/index.yaml

# Sync asset classes and recipes
refgenie data-channel sync my-channel --exists-ok

# List available asset classes
refgenie1 asset_class list

# Build an asset using your recipe
refgenie1 build genome/my_index --files fasta=/path/to/genome.fa
```

## Best Practices

1. **Version Control**: Always version your asset classes and recipes using semantic versioning
2. **Documentation**: Include clear descriptions in all YAML files
3. **Validation**: Run validation scripts before publishing updates
4. **Testing**: Test recipes locally before publishing
5. **Backwards Compatibility**: Avoid breaking changes to existing definitions
6. **Naming Conventions**: Use consistent, descriptive names for assets and recipes
7. **Docker Images**: Use specific version tags for Docker images, not `latest`
8. **Security**: Never include sensitive information in public channels

## Troubleshooting

### Common Issues

**Index file not updating:**

- Ensure `build_index.py` is run after adding/modifying files
- Check GitHub Actions logs if using CI/CD

**404 errors when accessing channel:**

- Verify GitHub Pages is enabled and deployed
- Check the exact URL structure matches your repository name
- Ensure all files are committed and pushed

**YAML validation errors:**

- Use a YAML linter to check syntax
- Ensure all required fields are present
- Check for proper indentation (spaces, not tabs)

**Recipe execution failures:**

- Verify Docker image is accessible
- Check command templates for syntax errors
- Ensure all template variables are properly escaped

## Contributing to the Official Channel

To contribute to the official refgenie recipes channel:

1. Fork the repository: https://github.com/refgenie/recipes
2. Add your asset classes and recipes
3. Run validation: `python data_channel_check.py .`
4. Submit a pull request with a clear description

The refgenie team will review and merge approved contributions.

## Additional Resources

- [Refgenie Data Channels Documentation](data_channels.md)
- [Refgenieserver Data Channels](server/data_channels.md)
- [Official Recipes Repository](https://github.com/refgenie/recipes)