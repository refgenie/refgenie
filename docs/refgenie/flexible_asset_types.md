# Asset classes and recipes

Refgenie describes each kind of asset with two YAML files, and anyone can write them.

- An **asset class** says what an asset *contains*: its seek keys (the files or values you can look up) and how a server offers it for download. Examples: `fasta`, `gtf`, `bwa_index`. See the [asset class specification](asset_class_specification.md).
- A **recipe** says how to *build* an asset of one class: its input assets, files, and parameters, the commands to run, and the container image that has the tools. See the [recipe specification](recipe_specification.md).

Classes and recipes are separate, so several recipes can build the same class. Recipe inputs ask for an asset *class*, not one particular asset, so a genome can hold several groups of one class (for example `fasta` and `fasta_txome`), and a build can use any of them as input. See [Asset groups and naming](naming_assets.md).

Asset classes and recipes are shared through [data channels](data_channels.md). The usual workflow:

1. Write an asset class for the new kind of asset.
2. Write a recipe that builds it.
3. Publish both in a data channel ([how](set_up_data_channel.md)).
4. Others sync the channel and run `refgenie build`.
