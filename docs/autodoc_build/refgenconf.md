<script>
document.addEventListener('DOMContentLoaded', (event) => {
  document.querySelectorAll('h3 code').forEach((block) => {
    hljs.highlightBlock(block);
  });
});
</script>

<style>
h3 { 
    padding-left: 22px;
    text-indent: -15px;
 }
h3 .hljs {
    padding-left: 20px;
    margin-left: 0px;
    text-indent: -15px;
    martin-bottom: 0px;
}
h4, table, p, li{ margin-left: 30px; }
h4 { 
    font-style: italic;
    font-size: 1em;
    margin-bottom: 0px;
}

</style>
# Package `refgenconf` Documentation

## <a name="RefGenConf"></a> Class `RefGenConf`
A sort of oracle of available reference genome assembly assets


```python
def assets_dict(self)
```

Map each assembly name to a list of available asset names.
#### Returns:

- `Mapping[str, Iterable[str]]`:  mapping from assembly name tocollection of available asset names.




```python
def assets_str(self, offset_text='  ', asset_sep='; ', genome_assets_delim=': ')
```

Create a block of text representing genome-to-asset mapping.
#### Parameters:

- `offset_text` (`str`):  text that begins each line of the textrepresentation that's produced
- `asset_sep` (`str`):  the delimiter between names of types of assets,within each genome line
- `genome_assets_delim` (`str`):  the delimiter to place betweenreference genome assembly name and its list of asset names


#### Returns:

- `str`:  text representing genome-to-asset mapping




```python
def genomes_list(self)
```

Get a list of this configuration's reference genome assembly IDs.
#### Returns:

- `Iterable[str]`:  list of this configuration's reference genomeassembly IDs




```python
def genomes_str(self)
```

Get as single string this configuration's reference genome assembly IDs.
#### Returns:

- `str`:  single string that lists this configuration's knownreference genome assembly IDs




```python
def get_asset(self, genome_name, asset_name, strict_exists=True, check_exist=<function RefGenConf.<lambda> at 0x7f69a6069c80>)
```

Get an asset for a particular assembly.
#### Parameters:

- `genome_name` (`str`):  name of a reference genome assembly of interest
- `asset_name` (`str`):  name of the particular asset to fetch
- `strict_exists` (`bool | NoneType`):  how to handle case in whichpath doesn't exist; True to raise IOError, False to raise RuntimeWarning, and None to do nothing at all
- `check_exist` (`function(callable) -> bool`):  how to check forasset/path existence


#### Returns:

- `str`:  path to the asset


#### Raises:

- `TypeError`:  if the existence check is not a one-arg function
- `refgenconf.MissingGenomeError`:  if the named assembly isn't knownto this configuration instance
- `refgenconf.MissingAssetError`:  if the names assembly is known tothis configuration instance, but the requested asset is unknown




```python
def list_assets_by_genome(self, genome=None)
```

List types/names of assets that are available for one--or all--genomes.
#### Parameters:

- `genome` (`str | NoneType`):  reference genome assembly ID, optional;if omitted, the full mapping from genome to asset names


#### Returns:

- `Iterable[str] | Mapping[str, Iterable[str]]`:  collection ofasset type names available for particular reference assembly if one is provided, else the full mapping between assembly ID and collection available asset type names




```python
def list_genomes_by_asset(self, asset=None)
```

List assemblies for which a particular asset is available.
#### Parameters:

- `asset` (`str | NoneType`):  name of type of asset of interest, optional


#### Returns:

- `Iterable[str] | Mapping[str, Iterable[str]]`:  collection ofassemblies for which the given asset is available; if asset argument is omitted, the full mapping from name of asset type to collection of assembly names for which the asset key is available will be returned.




```python
def list_remote(self, get_url=<function RefGenConf.<lambda> at 0x7f69a6069ea0>)
```

List genomes and assets available remotely.
#### Parameters:

- `get_url` (`function(refgenconf.RefGenConf) -> str`):  how to determineURL request, given RefGenConf instance


#### Returns:

- `str, str`:  text reps of remotely available genomes and assets




```python
def pull_asset(self, genome, assets, genome_config, unpack=True, get_json_url=<function RefGenConf.<lambda> at 0x7f69a6064048>, get_main_url=None)
```

Download and possibly unpack one or more assets for a given ref gen.
#### Parameters:

- `genome` (`str`):  name of a reference genome assembly of interest
- `assets` (`str`):  name(s) of particular asset(s) to fetch
- `genome_config` (`str`):  path to genome configuration file to update
- `unpack` (`bool`):  whether to unpack a tarball
- `get_json_url` (`function(str, str, str) -> str`):  how to build URL fromgenome server URL base, genome, and asset
- `get_main_url` (`function(str) -> str`):  how to get archive URL frommain URL


#### Returns:

- `Iterable[(str, str | NoneType)]`:  collection of pairs of assetname and folder name (key-value pair with which genome config file is updated) if pull succeeds, else asset key and a null value.


#### Raises:

- `TypeError`:  if the assets argument is neither string nor otherIterable
- `refgenconf.UnboundEnvironmentVariablesError`:  if genome folderpath contains any env. var. that's unbound




```python
def update_genomes(self, genome, asset=None, data=None)
```

Updates the genomes in RefGenConf object at any level. If a requested genome-asset mapping is missing, it will be created
#### Parameters:

- `genome` (`str`):  genome to be added/updated
- `asset` (`str`):  asset to be added/updated
- `data` (`Mapping`):  data to be added/updated


#### Returns:

- `RefGenConf`:  updated object




## <a name="MissingAssetError"></a> Class `MissingAssetError`
Error type for request of an unavailable genome asset.


## <a name="MissingGenomeError"></a> Class `MissingGenomeError`
Error type for request of unknown genome/assembly.


```python
def select_genome_config(filename, conf_env_vars=None)
```

Get path to genome configuration file.
#### Parameters:

- `filename` (`str`):  name/path of genome configuration file
- `conf_env_vars` (`Iterable[str]`):  names of environment variables toconsider; basically, a prioritized search list


#### Returns:

- `str`:  path to genome configuration file





*Version Information: `refgenconf` v0.1.0, generated by `lucidoc` v0.4.0*