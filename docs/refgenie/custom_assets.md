# Use custom assets

Refgenie records every asset you pull or build. You can also have it manage files you made yourself, so they are found the same way as every other asset.

## Build a custom asset

The best option is to script the build as a [recipe](recipe_specification.md), with an [asset class](asset_class_specification.md) for its output. Refgenie then builds the asset and records how it was made. To share the recipe, publish it in a [data channel](set_up_data_channel.md).

## Add an existing folder

If the asset cannot be scripted, add it with `refgenie add`. Say you made an hg38 annotation by hand. Put its folder in your genome folder, then give the registry path (with an asset name), the folder's path *relative to the genome folder*, and an asset class:

```console
refgenie add hg38/manual_anno:v1 --path manual_anno_folder --asset-class manual_anno
```

The asset class must already be registered (`refgenie asset-class add --source manual_anno.yaml`), and each of its path seek keys must match a file in the folder. Refgenie moves the folder into its own storage.

Now the asset works like any other: `refgenie list` shows it, `refgenie seek hg38/manual_anno` returns its path, and so does `r.asset.seek(...)` in Python.

Use `-k`/`--seek-keys name=value` to attach non-path seek key values, and `-d`/`--description` to describe the asset.
