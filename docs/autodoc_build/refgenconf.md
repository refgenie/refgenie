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
def __init__(self, filepath=None, entries=None, writable=False, wait_max=10)
```

Create the config instance by with a filepath or key-value pairs.
#### Parameters:

- `filepath` (`str`):  a path to the YAML file to read
- `entries` (`Iterable[(str, object)] | Mapping[str, object]`): config filepath or collection of key-value pairs
- `writable` (`bool`):  whether to create the object with write capabilities
- `wait_max` (`int`):  how long to wait for creating an object when the file that data will be read from is locked


#### Raises:

- `refgenconf.MissingConfigDataError`:  if a required configurationitem is missing
- `ValueError`:  if entries is given as a string and is not a file




```python
def assets_dict(self, genome=None, order=None, include_tags=False)
```

Map each assembly name to a list of available asset names.
#### Parameters:

- `order` (`function(str) -> object`):  how to key genome IDs for sort
- `genome` (`list[str] | str`):  genomes that the assets should be found for
- `include_tags` (`bool`):  whether asset tags should be included in the returned dict


#### Returns:

- `Mapping[str, Iterable[str]]`:  mapping from assembly name tocollection of available asset names.




```python
def assets_str(self, offset_text='  ', asset_sep=', ', genome_assets_delim='/ ', genome=None, order=None)
```

Create a block of text representing genome-to-asset mapping.
#### Parameters:

- `offset_text` (`str`):  text that begins each line of the textrepresentation that's produced
- `asset_sep` (`str`):  the delimiter between names of types of assets,within each genome line
- `genome_assets_delim` (`str`):  the delimiter to place betweenreference genome assembly name and its list of asset names
- `genome` (`list[str] | str`):  genomes that the assets should be found for
- `order` (``):  function(str) -> object how to key genome IDs and assetnames for sort


#### Returns:

- `str`:  text representing genome-to-asset mapping




```python
def chk_digest_update_child(self, genome, remote_asset_name, child_name, server_url)
```

Check local asset digest against the remote one and populate children of the asset with the provided asset:tag.

In case the local asset does not exist, the config is populated with the remote asset digest and children data
#### Parameters:

- `genome` (`str`):  name of the genome to check the asset digests for
- `remote_asset_name` (`str`): tag
- `child_name` (`str`):  name to be appended to the children of the parent
- `server_url` (`str`):  address of the server to query for the digests


#### Raises:

- `RefgenconfError`:  if the local digest does not match its remote counterpart




```python
def file_path(self)
```

Return the path to the config file or None if not set
#### Returns:

- `str | None`:  path to the file the object will would to




```python
def filepath(self, genome, asset, tag, ext='.tgz', dir=False)
```

Determine path to a particular asset for a particular genome.
#### Parameters:

- `genome` (`str`):  reference genome ID
- `asset` (`str`):  asset name
- `tag` (`str`):  tag name
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
def get_asset(self, genome_name, asset_name, tag_name=None, seek_key=None, strict_exists=True, check_exist=<function RefGenConf.<lambda> at 0x10bdfa950>, enclosing_dir=False)
```

Get an asset for a particular assembly.
#### Parameters:

- `genome_name` (`str`):  name of a reference genome assembly of interest
- `asset_name` (`str`):  name of the particular asset to fetch
- `tag_name` (`str`):  name of the particular asset tag to fetch
- `seek_key` (`str`):  name of the particular subasset to fetch
- `strict_exists` (`bool | NoneType`):  how to handle case in whichpath doesn't exist; True to raise IOError, False to raise RuntimeWarning, and None to do nothing at all
- `check_exist` (`function(callable) -> bool`):  how to check forasset/path existence
- `enclosing_dir` (`bool`):  whether a path to the entire enclosing directory should be returned, e.g.for a fasta asset that has 3 seek_keys pointing to 3 files in an asset dir, that asset dir is returned


#### Returns:

- `str`:  path to the asset


#### Raises:

- `TypeError`:  if the existence check is not a one-arg function
- `refgenconf.MissingGenomeError`:  if the named assembly isn't knownto this configuration instance
- `refgenconf.MissingAssetError`:  if the names assembly is known tothis configuration instance, but the requested asset is unknown




```python
def get_asset_digest(self, genome, asset, tag=None)
```

Returns the digest for the specified asset. The defined default tag will be used if not provided as an argument
#### Parameters:

- `genome` (`str`):  genome identifier
- `asset` (`str`):  asset identifier
- `tag` (`str`):  tag identifier


#### Returns:

- `str`:  asset digest for the tag




```python
def get_default_tag(self, genome, asset, use_existing=True)
```

Determine the asset tag to use as default. The one indicated by the 'default_tag' key in the asset section is returned. If no 'default_tag' key is found, by default the first listed tag is returned with a RuntimeWarning. This behavior can be turned off with use_existing=False
#### Parameters:

- `genome` (`str`):  name of a reference genome assembly of interest
- `asset` (`str`):  name of the particular asset of interest
- `use_existing` (`bool`):  whether the first tag in the config should be returned in case there is no defaulttag defined for an asset


#### Returns:

- `str`:  name of the tag to use as the default one




```python
def get_genome_attributes(self, genome)
```

Get the dictionary attributes, like checksum, contents, description. Does not return the assets.
#### Parameters:

- `genome` (`str`):  genome to get the attributes dict for


#### Returns:

- `Mapping[str, str]`:  available genome attributes




```python
def initialize_config_file(self, filepath=None)
```

Initialize genome configuration file on disk
#### Parameters:

- `filepath` (`str`):  a valid path where the configuration file should be initialized


#### Returns:

- `str`:  the filepath the file was initialized at


#### Raises:

- `OSError`:  in case the file could not be initialized due to insufficient permissions or pre-existence
- `TypeError`:  if no valid filepath cat be determined




```python
def is_asset_complete(self, genome, asset, tag)
```

Check whether all required tag attributes are defined in the RefGenConf object. This is the way we determine tag completeness.
#### Parameters:

- `genome` (`str`):  genome to be checked
- `asset` (`str`):  asset package to be checked
- `tag` (`str`):  tag to be checked


#### Returns:

- `bool`:  the decision




```python
def list_assets_by_genome(self, genome=None, order=None, include_tags=False)
```

List types/names of assets that are available for one--or all--genomes.
#### Parameters:

- `genome` (`str | NoneType`):  reference genome assembly ID, optional;if omitted, the full mapping from genome to asset names
- `order` (``):  function(str) -> object how to key genome IDs and assetnames for sort
- `include_tags` (`bool`):  whether asset tags should be included in the returned dict


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
def list_local(self, genome=None, order=None)
```

List locally available reference genome IDs and assets by ID.
#### Parameters:

- `genome` (`list[str] | str`):  genomes that the assets should be found for
- `order` (``):  function(str) -> object how to key genome IDs and assetnames for sort


#### Returns:

- `str, str`:  text reps of locally available genomes and assets




```python
def list_remote(self, genome=None, order=None, get_url=<function RefGenConf.<lambda> at 0x10bdfad08>)
```

List genomes and assets available remotely.
#### Parameters:

- `get_url` (`function(refgenconf.RefGenConf) -> str`):  how to determineURL request, given RefGenConf instance
- `genome` (`list[str] | str`):  genomes that the assets should be found for
- `order` (``):  function(str) -> object how to key genome IDs and assetnames for sort


#### Returns:

- `str, str`:  text reps of remotely available genomes and assets




```python
def pull_asset(self, genome, asset, tag, unpack=True, force=None, get_json_url=<function RefGenConf.<lambda> at 0x10bdfaf28>, build_signal_handler=<function _handle_sigint at 0x10bd65620>)
```

Download and possibly unpack one or more assets for a given ref gen.
#### Parameters:

- `genome` (`str`):  name of a reference genome assembly of interest
- `asset` (`str`):  name of particular asset to fetch
- `tag` (`str`):  name of particular tag to fetch
- `unpack` (`bool`):  whether to unpack a tarball
- `force` (`bool | NoneType`):  how to handle case in which asset pathalready exists; null for prompt (on a per-asset basis), False to effectively auto-reply No to the prompt to replace existing file, and True to auto-replay Yes for existing asset replacement.
- `get_json_url` (`function(str, str) -> str`):  how to build URL fromgenome server URL base, genome, and asset
- `build_signal_handler` (`function(str) -> function`):  how to createa signal handler to use during the download; the single argument to this function factory is the download filepath
- `update` (`bool`):  whether the object should be updated with downloaded archive data


#### Returns:

- `(list[str], dict, str)`:  a list of genome, asset, tag namesand a key-value pair with which genome config file should be updated if pull succeeds, else asset key and a null value


#### Raises:

- `refgenconf.UnboundEnvironmentVariablesError`:  if genome folderpath contains any env. var. that's unbound
- `refgenconf.RefGenConfError`:  if the object update is requested ina non-writable state




```python
def remove_asset_from_relatives(self, genome, asset, tag)
```

Remove any relationship links associated with the selected asset
#### Parameters:

- `genome` (`str`):  genome to be removed from its relatives' relatives list
- `asset` (`str`):  asset to be removed from its relatives' relatives list
- `tag` (`str`):  tag to be removed from its relatives' relatives list




```python
def remove_assets(self, genome, asset, tag=None)
```

Remove data associated with a specified genome:asset:tag combination. If no tags are specified, the entire asset is removed from the genome.

If no more tags are defined for the selected genome:asset after tag removal,
the parent asset will be removed as well
If no more assets are defined for the selected genome after asset removal,
the parent genome will be removed as well
#### Parameters:

- `genome` (`str`):  genome to be removed
- `asset` (`str`):  asset package to be removed
- `tag` (`str`):  tag to be removed


#### Returns:

- `RefGenConf`:  updated object


#### Raises:

- `TypeError`:  if genome argument type is not a list or str




```python
def set_default_pointer(self, genome, asset, tag, force=False)
```

Point to the selected tag by default
#### Parameters:

- `genome` (`str`):  name of a reference genome assembly of interest
- `asset` (`str`):  name of the particular asset of interest
- `tag` (`str`):  name of the particular asset tag to point to by default
- `force` (`bool`):  whether the default tag change should be forced (even if it exists)




```python
def tag_asset(self, genome, asset, tag, new_tag)
```

Retags the asset selected by the tag with the new_tag. Prompts if default already exists and overrides upon confirmation.

This method does not override the original asset entry in the RefGenConf object. It creates its copy and tags
it with the new_tag.
Additionally, if the retagged asset has any children their parent will be retagged as new_tag that was
introduced upon this method execution.
#### Parameters:

- `genome` (`str`):  name of a reference genome assembly of interest
- `asset` (`str`):  name of particular asset of interest
- `tag` (`str`):  name of the tag that identifies the asset of interest
- `new_tag` (`str`):  name of particular the new tag


#### Returns:

- `bool`:  a logical indicating whether the tagging was successful


#### Raises:

- `ValueError`:  when the original tag is not specified




```python
def update_assets(self, genome, asset=None, data=None)
```

Updates the genomes in RefGenConf object at any level. If a requested genome-asset mapping is missing, it will be created
#### Parameters:

- `genome` (`str`):  genome to be added/updated
- `asset` (`str`):  asset to be added/updated
- `data` (`Mapping`):  data to be added/updated


#### Returns:

- `RefGenConf`:  updated object




```python
def update_genome_servers(self, url, reset=False)
```

Update the list of genome_servers.

Use reset argument to overwrite the current list. Otherwise the current one will be appended to.
#### Parameters:

- `url` (`list[str] | str`):  url(s) to update the genome_servers list with
- `reset` (`bool`):  whether the current list should be overwritten




```python
def update_genomes(self, genome, data=None)
```

Updates the genomes in RefGenConf object at any level. If a requested genome is missing, it will be added
#### Parameters:

- `genome` (`str`):  genome to be added/updated
- `data` (`Mapping`):  data to be added/updated


#### Returns:

- `RefGenConf`:  updated object




```python
def update_relatives_assets(self, genome, asset, tag=None, data=None, children=False)
```

A convenience method which wraps the update assets and uses it to update the asset relatives of an asset.
#### Parameters:

- `genome` (`str`):  genome to be added/updated
- `asset` (`str`):  asset to be added/updated
- `tag` (`str`):  tag to be added/updated
- `data` (`list`):  asset parents to be added/updated
- `children` (`bool`):  a logical indicating whether the relationship to be added is 'children'


#### Returns:

- `RefGenConf`:  updated object




```python
def update_seek_keys(self, genome, asset, tag=None, keys=None)
```

A convenience method which wraps the update assets and uses it to update the seek keys for a tagged asset.
#### Parameters:

- `genome` (`str`):  genome to be added/updated
- `asset` (`str`):  asset to be added/updated
- `tag` (`str`):  tag to be added/updated
- `keys` (`Mapping`):  seek_keys to be added/updated


#### Returns:

- `RefGenConf`:  updated object




```python
def update_tags(self, genome, asset=None, tag=None, data=None)
```

Updates the genomes in RefGenConf object at any level. If a requested genome-asset-tag mapping is missing, it will be created
#### Parameters:

- `genome` (`str`):  genome to be added/updated
- `asset` (`str`):  asset to be added/updated
- `tag` (`str`):  tag to be added/updated
- `data` (`Mapping`):  data to be added/updated


#### Returns:

- `RefGenConf`:  updated object




```python
def writable(self)
```

Return writability flag or None if not set
#### Returns:

- `bool | None`:  whether the object is writable now




## <a name="GenomeConfigFormatError"></a> Class `GenomeConfigFormatError`
Exception for invalid genome config file format.


```python
def __init__(self, msg)
```

Initialize self.  See help(type(self)) for accurate signature.



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
def select_genome_config(filename=None, conf_env_vars=['REFGENIE'], **kwargs)
```

Get path to genome configuration file.
#### Parameters:

- `filename` (`str`):  name/path of genome configuration file
- `conf_env_vars` (`Iterable[str]`):  names of environment variables toconsider; basically, a prioritized search list


#### Returns:

- `str`:  path to genome configuration file







*Version Information: `refgenconf` v0.6.1, generated by `lucidoc` v0.4.1*