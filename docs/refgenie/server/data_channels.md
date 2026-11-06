# Serve data channels

A [data channel](../data_channels.md) is a set of asset class and recipe files listed in an `index.yaml`, hosted as static files. The refgenie server re-serves the channels in its config file, so a client can get recipes and asset classes from the server it already uses. Each channel keeps its own name and layout: the server returns the same `index.yaml` and the same relative paths as the upstream channel.

## Configure channels

List the channels in a YAML file:

```yaml
registry:
  name: registry
  protocol: https
  index_url: https://refgenie.github.io/refgenie-registry/index.yaml

my_lab:
  name: my_lab
  protocol: https
  index_url: https://example.org/my-recipes/index.yaml
```

- `name`: the channel name, used in URLs. It cannot contain `__`.
- `protocol`: `http` or `https`. It must match the scheme of `index_url`.
- `index_url`: the URL of the channel's `index.yaml`.

Point the server at the file before starting it:

```bash
export DATA_CHANNELS_CONFIG_PATH=/path/to/data_channels.yaml
refgenie serve
```

Without `DATA_CHANNELS_CONFIG_PATH`, the server uses the packaged `refgenie/server/data_channels.yaml`, which lists the `registry` channel. The server will not start if the file is missing, empty or invalid.

## Endpoints

| Request | Response |
|---|---|
| `GET /v4/data_channel/` | The configured channels: name, protocol and `index_url` |
| `GET /v4/data_channel/{channel}/index.yaml` | That channel's index, fetched from upstream |
| `GET /v4/data_channel/{channel}/{path}` | A redirect to the file at `path` on the upstream channel |

`path` is a relative path from the channel's `index.yaml`, such as `recipes/bwa_index/recipe.yaml`. The server redirects only to files inside the upstream index's own folder; other paths get a 404.

To publish a channel of your own, see [Set up your own data channel](../set_up_data_channel.md).
