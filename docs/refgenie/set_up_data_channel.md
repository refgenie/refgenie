# Set up your own data channel

A [data channel](data_channels.md) is a folder of asset class and recipe YAML files, plus an `index.yaml` that lists them, served at a URL. Publish one to share your own kinds of assets and build recipes with your group or with everyone.

## 1. Lay out the folder

```
my-data-channel/
    index.yaml                  # lists every file below
    asset_classes/
        my_index.yaml
    recipes/
        my_index/recipe.yaml
```

The [refgenie-registry](https://github.com/refgenie/refgenie-registry) repository uses this layout and is a good model. Its [`tools/build_channel.py`](https://github.com/refgenie/refgenie-registry/blob/master/tools/build_channel.py) writes `index.yaml` for you, and [`tools/validate_channel.py`](https://github.com/refgenie/refgenie-registry/blob/master/tools/validate_channel.py) checks a channel folder or URL.

## 2. Write an asset class

An asset class lists the files (seek keys) that make up a kind of asset, and how a server offers them for download.

```yaml
name: my_index
version: 0.0.1
description: Custom index format for my tool
serving_modes: [archive]
seek_keys:
  index_file:
    value: "{genome}.idx"
    description: Main index file
    type: file
  metadata:
    value: "{genome}.meta"
    description: Index metadata
    type: file
```

See the [asset class specification](asset_class_specification.md).

## 3. Write a recipe

A recipe says how to build an asset of that class.

```yaml
name: my_index
version: 0.0.1
output_asset_class: my_index
description: Build a custom index from a FASTA file
input_files: null
input_params: null
input_assets:
  fasta:
    asset_class: fasta
    description: fasta asset for genome
    default: fasta
docker_image: docker.io/myorg/my-indexer:1.2.0
command_templates:
  - my-indexer build {{values.genome_folder}}/{{values.assets["fasta"].seek_keys_dict["fasta"]}} -o {{values.output_folder}}/{{values.genome_digest}}.idx
  - my-indexer meta {{values.output_folder}}/{{values.genome_digest}}.idx > {{values.output_folder}}/{{values.genome_digest}}.meta
custom_seek_keys:
  version: "my-indexer --version | head -1"
default_asset: "{{values.custom_seek_keys.version}}"
```

See the [recipe specification](recipe_specification.md).

## 4. Write the index

```yaml
asset_class:
  dir: asset_classes
  files:
    - my_index.yaml
recipe:
  dir: recipes
  files:
    - my_index/recipe.yaml
```

Each file path is relative to its section's `dir`, which is relative to the folder that holds `index.yaml`.

## 5. Host it

Any static web host works. With GitHub Pages, push the folder to a repository, turn on Pages in the repository settings, and the channel is at `https://<username>.github.io/<repository>/index.yaml`. The refgenie-registry repository has a [workflow](https://github.com/refgenie/refgenie-registry/blob/master/.github/workflows/publish-channel.yaml) that rebuilds and publishes its channel on every push.

You can also use a channel straight from disk, with the `local` type and the path to `index.yaml`.

Instead of static files, a server can answer the same requests: `GET <base>/index.yaml` returns the index, and `GET <base>/<dir>/<file>` returns each YAML file. A [refgenie server](server/data_channels.md) can re-serve channels this way.

## 6. Test it

```bash
refgenie data-channel add my-channel https https://<username>.github.io/<repository>/index.yaml
refgenie data-channel validate my-channel
refgenie data-channel sync my-channel --exists-ok
refgenie asset-class list
refgenie build hg38/my_index
```

## Tips

- Use [semantic versions](https://semver.org/) for asset classes and recipes; refgenie rejects other version strings. Give a changed file a new version instead of editing an old one.
- Pin container images to a version tag, not `latest`.
- Test each recipe with `refgenie build --container` before you publish it.
- To add to the refgenie channel, open a pull request on [refgenie-registry](https://github.com/refgenie/refgenie-registry); its [contributing guide](https://github.com/refgenie/refgenie-registry/blob/master/CONTRIBUTING.md) explains the checks.

## See also

- [Data channels](data_channels.md)
- [Serve data channels](server/data_channels.md)
