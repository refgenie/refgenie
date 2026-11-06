# FAQ

## Why isn't the refgenie command on my PATH?

`pip install --user` puts Python command-line programs in `~/.local/bin`. Add that folder to your path with `export PATH=$PATH:~/.local/bin`, and put the line in your `.bashrc` to keep it.

## Can I use refgenie with genome files I already have?

Yes. Register them with [refgenie add](custom_assets.md), and then use `seek`, the Python API, and the rest of refgenie with them.

## Can I add an asset to the public server?

You can run your own server with `refgenie serve` and host any assets you like; see [Run a server](server_tutorial.md). To add a new kind of asset for everyone, publish its recipe in a [data channel](set_up_data_channel.md). If you think an asset belongs on the public server, [open an issue](https://github.com/refgenie/refgenie/issues).

## Can I keep several versions of an asset?

Yes. Each asset group (such as `bowtie2_index`) can hold several named assets, such as `2.5.5` and `2.4.0`, and one of them is the default. See [Managing asset groups and naming](naming_assets.md).

## Can several people share one refgenie setup?

Yes. Point everyone at the same database and genome folder. For several users writing at once, use PostgreSQL; see [Database backends](database.md).

## How can I tell how an asset was built?

Each asset records how it was built: the recipe and its version, the input assets and files, the parameters, the container image, and the refgenie version. A server returns this record with each asset from its `/v4/assets` endpoint.
