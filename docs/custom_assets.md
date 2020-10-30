# Use custom assets

Refgenie will write the genome configuration file automatically for any assets that you obtain using either `refgenie pull` or `refgenie build` -- but you can also take advantage of the refgenie system to manage your own custom assets in a few ways.

## Build a custom asset

The preferred option would be to script your asset building and then allow `refgenie` to manage it. In the next major version of `refgenie`, we intend to allow the `build` function to build arbitrary assets. So, all you would need to do is be able to provide a scripted recipe and you could use `refgenie` to build and manage those assets automatically. In the meantime, or if you have assets that you want managed but *can't* be scripted for some reason...

## Add custom assets

You can add additional assets with the `refgenie add` command. Just provide the genome, asset, and path *relative to the genome folder*. What this means is that you can exploit the refgenie system to manage and access your own assets. For example, say you have an hg38 annotation called *manual_annotation*, which you produced by hand. You can simply put that in your genomes folder (under `hg38/annotation_folder_dir`), and then add an entry to your genome configuration file:

```console
refgenie add hg38/manual_anno --path annotation_folder_dir
```

If you want to, you could also just edit the config file by hand by adding this kind of information:

```yaml
genomes:
  511fb1178275e7d529560d53b949dba40815f195623bce8e:
    aliases:
      - hg38
      - human
    assets:
      manual_anno:
        tags:
          default:
            asset_path: manual_anno
            asset_description: Manual annotations from project X
            seek_keys:
              anno1: anno1.txt
              anno2: anno2.txt
        default_tag: default
```

The refgenie-compatible genome digest can be determined this way:

```python
from refgenconf.seqcol import SeqColClient
digest, _ = SeqColClient({}).load_fasta("path/hg38.fa")
# digest -> 511fb1178275e7d529560d53b949dba40815f195623bce8e
```

Now, you can access this asset with `refgenie` the same way you do all other assets... `refgenie list` will include it, `refgenie seek -g gh38 -a manual_anno` will retrieve the path, and from within python, `RefGenConf.seek('hg38', 'manual_anno')` will also work. The advantage of doing this is that it lets you include *all* your genome-associated resources, including manual ones, within the same framework.
