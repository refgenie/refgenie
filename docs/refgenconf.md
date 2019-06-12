# Refgenie from within Python

One of refgenie's key strengths is that independent Python tools can rely on our Python object that provides an API for use by Python tool developers. For this we have a separate Python package called `refgenconf` which provides a data type with which to access local and remote genome assets.

## Installing

You should already have `refgenconf` if you've installed `refgenie`, but if needed you can also install it separately with some variant of `pip install refgenconf`, depending on your usage context.

## Quick start
Create a `RefGenConf` object, which is the package's main data type. You just need to give it a genome configuration file (in YAML format).

```python
import refgenconf
rgc = refgenconf.RefGenConf("genome_config.yaml")
```

Now, you can interact with it:
```python
print(rgc)
```

Use this to show all available remote assets:
```python
rgc.list_remote()
```

See the complete [refgenie python API](/autodoc_build/refgenconf) for details.


