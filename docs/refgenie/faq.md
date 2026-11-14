# FAQ

## Why isn't the refgenie command on my PATH?

`pip install --user` puts Python command-line programs in `~/.local/bin`. Add that folder to your path with `export PATH=$PATH:~/.local/bin`, and put the line in your `.bashrc` to keep it.

## Can I use refgenie with genome files I already have?

Yes. Register them with [refgenie add](custom_assets.md), and then use `seek`, the Python API, and the rest of refgenie with them.

## How do I ask for a genome or asset refgenie does not have?

- A genome refgenie does not list: [open a genome proposal](https://github.com/refgenie/refgenie-registry/issues/new?template=genome_proposal.yml).
- A listed genome missing an asset: [open an asset request](https://github.com/refgenie/refgenie-registry/issues/new?template=asset_request.yml). Put the genome's alias and digest in the form.
- Need it now: build it yourself; see [Build assets](build.md). To host assets of your own, run your own server with `refgenie serve` ([Run a server](server_tutorial.md)), and publish new recipes in a [data channel](set_up_data_channel.md).

## Where do I get the aligner for an index?

Refgenie gives you the index, not the aligner. The asset name is the tool version that built it: `refgenie seek rCRSd/bowtie2_index.version` prints `2.5.5`, so install bowtie2 2.5.5 and match the version. Use conda (`conda install -c bioconda bowtie2=2.5.5`), or, without conda, the prebuilt binary from the tool's GitHub releases, such as the [bowtie2 releases](https://github.com/BenLangmead/bowtie2/releases).

## Can I keep several versions of an asset?

Yes. Each asset group (such as `bowtie2_index`) can hold several named assets, such as `2.5.5` and `2.4.0`, and one of them is the default. See [Managing asset groups and naming](naming_assets.md).

## Can several people share one refgenie setup?

Yes. Point everyone at the same database and genome folder. For several users writing at once, use PostgreSQL; see [Database backends](database.md).

## How can I tell how an asset was built?

Each asset records how it was built: the recipe and its version, the input assets and files, the parameters, the container image, and the refgenie version. A server returns this record with each asset from its `/v4/assets` endpoint.
