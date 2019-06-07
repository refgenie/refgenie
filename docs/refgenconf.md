# Refgenie from within python

One of refgenie's key strengths is that independent python tools can rely on our python object that provides an API for use by python tool developers. For this we have a separate python package called `refgenconf` which provides a configuration object giving access to local and remote genome assets.

## Installing

You should already have `refgenconf` if you've installed `refgenie`, but you can also install it separately with `pip install refgenconf` if necessary.

## Quick start

Create a `RefGenConf` object, which is the primary object type from the python API. You just need to give it a yaml genome configuration file:

```
import refgenconf
rgc = refgenconf.RefGenConf("genome_config.yaml")
```

Now, you can interact with it:

```
print(rgc)
```

Use this to show all available remote assets:

```
rgc.list_remote()
```

See the complete [refgenie python API](/autodoc_build/refgenconf) for details.


