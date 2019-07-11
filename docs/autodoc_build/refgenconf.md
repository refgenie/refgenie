<script>
document.addEventListener('DOMContentLoaded', (event) => {
  document.querySelectorAll('h3 code').forEach((block) => {
    hljs.highlightBlock(block);
  });
});
</script>

<style>
h3 .content { 
    padding-left: 22px;
    text-indent: -15px;
 }
h3 .hljs .content {
    padding-left: 20px;
    margin-left: 0px;
    text-indent: -15px;
    martin-bottom: 0px;
}
h4 .content, table .content, p .content, li .content { margin-left: 30px; }
h4 .content { 
    font-style: italic;
    font-size: 1em;
    margin-bottom: 0px;
}

</style>


# Package `refgenconf` Documentation

## <a name="RefGenConf"></a> Class `RefGenConf`
A sort of oracle of available reference genome assembly assets


```python
def assets_dict(self, order=None)
```

Map each assembly name to a list of available asset names.
#### Parameters:

- `order` (``):  function(str) -> object how to key genome IDs for sort


#### Returns:

- `Mapping[str, Iterable[str]]`:  mapping from assembly name tocollection of available asset names.




```python
def assets_str(self, offset_text='  ', asset_sep='; ', genome_assets_delim=': ', order=None)
```

Create a block of text representing genome-to-asset mapping.
#### Parameters:

- `offset_text` (`str`):  text that begins each line of the textrepresentation that's produced
- `asset_sep` (`str`):  the delimiter between names of types of assets,within each genome line
- `genome_assets_delim` (`str`):  the delimiter to place betweenreference genome assembly name and its list of asset names
- `order` (``):  function(str) -> object how to key genome IDs and assetnames for sort


#### Returns:

- `str`:  text representing genome-to-asset mapping




```python
def filepath(self, genome, asset, ext='.tar')
```

Determine path to a particular asset for a particular genome.
#### Parameters:

- `genome` (`str`):  reference genome iD
- `asset` (`str`):  asset name
- `ext` (`str`):  file extension


#### Returns:

- `str`:  path to asset for given genome and asset kind/name




```python
def genomes_list(self, order=None)
```

Get a list of this configuration's reference genome assembly IDs.
#### Returns:

- `Iterable[str]`:  list of this configuration's reference genomeassembly IDs




```python
def genomes_str(self, order=None)
```

Get as single string this configuration's reference genome assembly IDs.
#### Parameters:

- `order` (``):  function(str) -> object how to key genome IDs for sort


#### Returns:

- `str`:  single string that lists this configuration's knownreference genome assembly IDs




```python
def get_asset(self, genome_name, asset_name, strict_exists=True, check_exist=<function RefGenConf.<lambda> at 0x7ff741bea268>)
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
def list_assets_by_genome(self, genome=None, order=None)
```

List types/names of assets that are available for one--or all--genomes.
#### Parameters:

- `genome` (`str | NoneType`):  reference genome assembly ID, optional;if omitted, the full mapping from genome to asset names
- `order` (``):  function(str) -> object how to key genome IDs and assetnames for sort


#### Returns:

- `Iterable[str] | Mapping[str, Iterable[str]]`:  collection ofasset type names available for particular reference assembly if one is provided, else the full mapping between assembly ID and collection available asset type names




```python
def list_genomes_by_asset(self, asset=None, order=None)
```

List assemblies for which a particular asset is available.
#### Parameters:

- `asset` (`str | NoneType`):  name of type of asset of interest, optional
- `order` (``):  function(str) -> object how to key genome IDs and assetnames for sort


#### Returns:

- `Iterable[str] | Mapping[str, Iterable[str]]`:  collection ofassemblies for which the given asset is available; if asset argument is omitted, the full mapping from name of asset type to collection of assembly names for which the asset key is available will be returned.




```python
def list_local(self, order=None)
```

List locally available reference genome IDs and assets by ID.
#### Parameters:

- `order` (``):  function(str) -> object how to key genome IDs and assetnames for sort


#### Returns:

- `str, str`:  text reps of locally available genomes and assets




```python
def list_remote(self, get_url=<function RefGenConf.<lambda> at 0x7ff741bea510>, order=None)
```

List genomes and assets available remotely.
#### Parameters:

- `get_url` (`function(refgenconf.RefGenConf) -> str`):  how to determineURL request, given RefGenConf instance
- `order` (``):  function(str) -> object how to key genome IDs and assetnames for sort


#### Returns:

- `str, str`:  text reps of remotely available genomes and assets




```python
def pull_asset(self, genome, assets, genome_config, unpack=True, force=None, get_json_url=<function RefGenConf.<lambda> at 0x7ff741bea620>, get_main_url=None, build_signal_handler=<function _handle_sigint at 0x7ff741c6ba60>)
```

Download and possibly unpack one or more assets for a given ref gen.
#### Parameters:

- `genome` (`str`):  name of a reference genome assembly of interest
- `assets` (`str`):  name(s) of particular asset(s) to fetch
- `genome_config` (`str`):  path to genome configuration file to update
- `unpack` (`bool`):  whether to unpack a tarball
- `force` (`bool | NoneType`):  how to handle case in which asset pathalready exists; null for prompt (on a per-asset basis), False to effectively auto-reply No to the prompt to replace existing file, and True to auto-replay Yes for existing asset replacement.
- `get_json_url` (`function(str, str, str) -> str`):  how to build URL fromgenome server URL base, genome, and asset
- `get_main_url` (`function(str) -> str`):  how to get archive URL frommain URL
- `build_signal_handler` (`function(str) -> function`):  how to createa signal handler to use during the download; the single argument to this function factory is the download filepath


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




## <a name="GenomeConfigFormatError"></a> Class `GenomeConfigFormatError`
Exception for invalid genome config file format.


## <a name="MissingAssetError"></a> Class `MissingAssetError`
Error type for request of an unavailable genome asset.


## <a name="MissingConfigDataError"></a> Class `MissingConfigDataError`
Missing required configuration instance items


## <a name="MissingGenomeError"></a> Class `MissingGenomeError`
Error type for request of unknown genome/assembly.


## <a name="RefgenconfError"></a> Class `RefgenconfError`
Base exception type for this package


## <a name="UnboundEnvironmentVariablesError"></a> Class `UnboundEnvironmentVariablesError`
Use of environment variable that isn't bound to a value.


```python
def select_genome_config(filename, conf_env_vars=None, **kwargs)
```

Get path to genome configuration file.
#### Parameters:

- `filename` (`str`):  name/path of genome configuration file
- `conf_env_vars` (`Iterable[str]`):  names of environment variables toconsider; basically, a prioritized search list


#### Returns:

- `str`:  path to genome configuration file







*Version Information: `refgenconf` v0.2.0, generated by `lucidoc` v0.4.0*