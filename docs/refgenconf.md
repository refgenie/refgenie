# Refgenie from within Python

Third-party python tools can rely on our Python object for access to refgenie assets. For this we have a Python package called `refgenconf` which provides a class with methods to access local and remote genome assets.

## Installing

You should already have `refgenconf` if you've installed `refgenie`, but if needed you can also install it separately with some variant of `pip install refgenconf`.

## Quick start
Create a `RefGenConf` object, which is the package's main data type. You just need to give it a refgenie genome configuration file (in YAML format). You can create a template using `refgenie init`.

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

In a tool, you're probably most interested in using refgenie to locate reference genome assets, for which you want to use the `get_asset` function. For example:

```python
# identify genome (perhaps provided by user)
genome = "hg38"

# get the local path to bowtie2 indexes:
bt2idx = rgc.get_asset(genome, "bowtie2_index")

# run bowtie2...
```

This enables you to write python software that will work on any computing environment without having to worry about passing around brittle environment-specific file paths.

See the complete [refgenconf python API](/autodoc_build/refgenconf) for more details.
