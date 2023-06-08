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
def __init__(self, filepath=None, entries=None, writable=False, wait_max=60, skip_read_lock=False, genome_exact=False, schema_source=None)
```

Create the config instance by with a filepath or key-value pairs.
#### Parameters:

- `filepath` (`str`):  a path to the YAML file to read
- `entries` (`Iterable[(str, object)] | Mapping[str, object]`): config filepath or collection of key-value pairs
- `writable` (`bool`):  whether to create the object with write capabilities
- `wait_max` (`int`):  how long to wait for creating an object when thefile that data will be read from is locked
- `skip_read_lock` (`bool`):  whether the file should not be locked forreading when object is created in read only mode


#### Raises:

- `refgenconf.MissingConfigDataError`:  if a required configurationitem is missing
- `ValueError`:  if entries is given as a string and is not a file




```python
def add(self, path, genome, asset, tag=None, seek_keys=None, force=False)
```

Add an external asset to the config
#### Parameters:

- `path` (`str`):  a path to the asset to add; must exist and be relativeto the genome_folder
- `genome` (`str`):  genome name
- `asset` (`str`):  asset name
- `tag` (`str`):  tag name
- `seek_keys` (`dict`):  seek keys to add
- `force` (`bool`):  whether to force existing asset overwrite




```python
def alias_dir(self)
```

Path to the genome alias directory
#### Returns:

- `str`:  path to the directory where the assets are stored




```python
def assets_str(self, offset_text='  ', asset_sep=', ', genome_assets_delim='/ ', genome=None, order=None)
```

Create a block of text representing genome-to-asset mapping.
#### Parameters:

- `offset_text` (`str`):  text that begins each line of the textrepresentation that's produced
- `asset_sep` (`str`):  the delimiter between names of types of assets,within each genome line
- `genome_assets_delim` (`str`):  the delimiter to place betweenreference genome assembly name and its list of asset names
- `genome` (`list[str] | str`):  genomes that the assets should be found for
- `order` (`function(str) -> object`):  how to key genome IDs and assetnames for sort


#### Returns:

- `str`:  text representing genome-to-asset mapping




```python
def cfg_remove_assets(self, genome, asset, tag=None, relationships=True)
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
- `relationships` (`bool`):  whether the asset being removed shouldbe removed from its relatives as well


#### Returns:

- `RefGenConf`:  updated object


#### Raises:

- `TypeError`:  if genome argument type is not a list or str




```python
def cfg_tag_asset(self, genome, asset, tag, new_tag, force=False)
```

Retags the asset selected by the tag with the new_tag. Prompts if default already exists and overrides upon confirmation.

This method does not override the original asset entry in the
RefGenConf object. It creates its copy and tags it with the new_tag.
Additionally, if the retagged asset has any children their parent will
 be retagged as new_tag that was introduced upon this method execution.
#### Parameters:

- `genome` (`str`):  name of a reference genome assembly of interest
- `asset` (`str`):  name of particular asset of interest
- `tag` (`str`):  name of the tag that identifies the asset of interest
- `new_tag` (`str`):  name of particular the new tag
- `force` (`bool`):  force any actions that require approval


#### Returns:

- `bool`:  a logical indicating whether the tagging was successful


#### Raises:

- `ValueError`:  when the original tag is not specified




```python
def chk_digest_update_child(self, genome, remote_asset_name, child_name, server_url)
```

Check local asset digest against the remote one and populate children of the asset with the provided asset:tag.

In case the local asset does not exist, the config is populated with the remote
 asset digest and children data
#### Parameters:

- `genome` (`str`):  name of the genome to check the asset digests for
- `remote_asset_name` (`str`): tag
- `child_name` (`str`):  name to be appended to the children of the parent
- `server_url` (`str`):  address of the server to query for the digests


#### Raises:

- `RefgenconfError`:  if the local digest does not match its remote counterpart




```python
def compare(self, genome1, genome2, explain=False)
```

Check genomes compatibility level. Compares Annotated Sequence Digests (ASDs) -- digested sequences and metadata
#### Parameters:

- `genome1` (`str`):  name of the first genome to compare
- `genome2` (`str`):  name of the first genome to compare
- `explain` (`bool`):  whether the returned code explanation shouldbe displayed


#### Returns:

- `int`:  compatibility code




```python
def data_dir(self)
```

Path to the genome data directory
#### Returns:

- `str`:  path to the directory where the assets are stored




```python
def file_path(self)
```

Path to the genome configuration file
#### Returns:

- `str`:  path to the genome configuration file




```python
def filepath(self, genome, asset, tag, ext='.tgz', dir=False)
```

Determine path to a particular asset for a particular genome.
#### Parameters:

- `genome` (`str`):  reference genome ID
- `asset` (`str`):  asset name
- `tag` (`str`):  tag name
- `ext` (`str`):  file extension
- `dir` (`bool`):  whether to return the enclosing directory instead of the file


#### Returns:

- `str`:  path to asset for given genome and asset kind/name




```python
def genome_aliases(self)
```

Mapping of human-readable genome identifiers to genome identifiers
#### Returns:

- `dict`:  mapping of human-readable genome identifiers to genomeidentifiers




```python
def genome_aliases_table(self)
```

Mapping of human-readable genome identifiers to genome identifiers
#### Returns:

- `dict`:  mapping of human-readable genome identifiers to genomeidentifiers




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

- `order` (`function(str) -> object`):  how to key genome IDs for sort


#### Returns:

- `str`:  single string that lists this configuration's knownreference genome assembly IDs




```python
def get_asds_path(self, genome)
```

Get path to the Annotated Sequence Digests JSON file for a given genome. Note that the path and/or genome may not exist.
#### Parameters:

- `genome` (`str`):  genome name


#### Returns:

- `str`:  ASDs path




```python
def get_asset_table(self, genomes=None, server_url=None, get_json_url=<function RefGenConf.<lambda> at 0x7f810efbc940>)
```

Get a rich.Table object representing assets available locally
#### Parameters:

- `genomes` (`list[str]`):  genomes to restrict the results with
- `server_url` (`str`):  server URL to query for the remote genome data
- `get_json_url` (`function(str, str) -> str`):  how to build URL fromgenome server URL base, genome, and asset


#### Returns:

- `rich.table.Table`:  table of assets available locally




```python
def get_default_tag(self, genome, asset, use_existing=True)
```

Determine the asset tag to use as default. The one indicated by the 'default_tag' key in the asset section is returned. If no 'default_tag' key is found, by default the first listed tag is returned with a RuntimeWarning. This behavior can be turned off with use_existing=False
#### Parameters:

- `genome` (`str`):  name of a reference genome assembly of interest
- `asset` (`str`):  name of the particular asset of interest
- `use_existing` (`bool`):  whether the first tag in the config should bereturned in case there is no default tag defined for an asset


#### Returns:

- `str`:  name of the tag to use as the default one




```python
def get_genome_alias(self, digest, fallback=False, all_aliases=False)
```

Get the human readable alias for a genome digest
#### Parameters:

- `digest` (`str`):  digest to find human-readable alias for
- `fallback` (`bool`):  whether to return the query digest in caseof failure
- `all_aliases` (`bool`):  whether to return all aliases instead of justthe first one


#### Returns:

- `str | list[str]`:  human-readable aliases


#### Raises:

- `GenomeConfigFormatError`:  if "genome_digests" section doesnot exist in the config
- `UndefinedAliasError`:  if a no alias has been defined for therequested digest




```python
def get_genome_alias_digest(self, alias, fallback=False)
```

Get the human readable alias for a genome digest
#### Parameters:

- `alias` (`str`):  alias to find digest for
- `fallback` (`bool`):  whether to return the query alias in caseof failure and in case it is one of the digests


#### Returns:

- `str`:  genome digest


#### Raises:

- `UndefinedAliasError`:  if the specified alias has been assigned toany digests




```python
def get_genome_attributes(self, genome)
```

Get the dictionary attributes, like checksum, contents, description. Does not return the assets.
#### Parameters:

- `genome` (`str`):  genome to get the attributes dict for


#### Returns:

- `Mapping[str, str]`:  available genome attributes




```python
def get_local_data_str(self, genome=None, order=None)
```

List locally available reference genome IDs and assets by ID.
#### Parameters:

- `genome` (`list[str] | str`):  genomes that the assets should be found for
- `order` (`function(str) -> object`):  how to key genome IDs and assetnames for sort


#### Returns:

- `str, str`:  text reps of locally available genomes and assets




```python
def get_remote_data_str(self, genome=None, order=None, get_url=<function RefGenConf.<lambda> at 0x7f810efb7670>)
```

List genomes and assets available remotely.
#### Parameters:

- `get_url` (`function(serverUrl, operationId) -> str`):  how to determineURL request, given server URL and endpoint operationID
- `genome` (`list[str] | str`):  genomes that the assets should be found for
- `order` (`function(str) -> object`):  how to key genome IDs and assetnames for sort


#### Returns:

- `str, str`:  text reps of remotely available genomes and assets




```python
def get_symlink_paths(self, genome, asset=None, tag=None, all_aliases=False)
```

Get path to the alias directory for the selected genome-asset-tag
#### Parameters:

- `genome` (`str`):  reference genome ID
- `asset` (`str`):  asset name
- `tag` (`str`):  tag name
- `all_aliases` (`bool`):  whether to return a collection of symboliclinks or just the first one from the alias list


#### Returns:

- `dict`: 




```python
def getseq(self, genome, locus, as_str=False)
```

Return the sequence found in a selected range and chromosome. Something like the refget protocol.
#### Parameters:

- `genome` (`str`):  name of the sequence identifier
- `locus` (`str`): 1-10'
- `as_str` (`bool`):  whether to convert the resurned object to stringand return just the sequence


#### Returns:

- `str | pyfaidx.FastaRecord | pyfaidx.Sequence`:  selected sequence




```python
def id(self, genome, asset, tag=None)
```

Returns the digest for the specified asset. The defined default tag will be used if not provided as an argument
#### Parameters:

- `genome` (`str`):  genome identifier
- `asset` (`str`):  asset identifier
- `tag` (`str`):  tag identifier


#### Returns:

- `str`:  asset digest for the tag




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
def initialize_genome(self, fasta_path, alias, fasta_unzipped=False, skip_alias_write=False)
```

Initialize a genome

Create a JSON file with Annotated Sequence Digests (ASDs)
for the FASTA file in the genome directory.
#### Parameters:

- `fasta_path` (`str`):  path to a FASTA file to initialize genome with
- `alias` (`str`):  alias to set for the genome
- `skip_alias_write` (`bool`):  whether to skip writing the alias to the file


#### Returns:

- `str, list[dict[]]`:  human-readable name for the genome




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
def list(self, genome=None, order=None, include_tags=False)
```

List local assets; map each namespace to a list of available asset names
#### Parameters:

- `order` (`callable(str) -> object`):  how to key genome IDs for sort
- `genome` (`list[str] | str`):  genomes that the assets should be found for
- `include_tags` (`bool`):  whether asset tags should be included in the returned dict


#### Returns:

- `Mapping[str, Iterable[str]]`:  mapping from assembly name tocollection of available asset names.




```python
def list_assets_by_genome(self, genome=None, order=None, include_tags=False)
```

List types/names of assets that are available for one--or all--genomes.
#### Parameters:

- `genome` (`str | NoneType`):  reference genome assembly ID, optional;if omitted, the full mapping from genome to asset names
- `order` (`function(str) -> object`):  how to key genome IDs and assetnames for sort
- `include_tags` (`bool`):  whether asset tags should be included in thereturned dict


#### Returns:

- `Iterable[str] | Mapping[str, Iterable[str]]`:  collection ofasset type names available for particular reference assembly if one is provided, else the full mapping between assembly ID and collection available asset type names




```python
def list_genomes_by_asset(self, asset=None, order=None)
```

List assemblies for which a particular asset is available.
#### Parameters:

- `asset` (`str | NoneType`):  name of type of asset of interest, optional
- `order` (`function(str) -> object`):  how to key genome IDs and assetnames for sort


#### Returns:

- `Iterable[str] | Mapping[str, Iterable[str]]`:  collection ofassemblies for which the given asset is available; if asset argument is omitted, the full mapping from name of asset type to collection of assembly names for which the asset key is available will be returned.




```python
def list_seek_keys_values(self, genomes=None, assets=None)
```

List values for all seek keys for the specified genome and asset. Leave the arguments out to get all seek keys values managed by refgenie.
#### Parameters:

- `genome_names` (`str | List[str]`):  optional list of genomes to include
- `asset_names` (`str | List[str]`):  optional list of assets to include


#### Returns:

- `dict`:  a nested dictionary with the seek key values




```python
def listr(self, genome=None, get_url=<function RefGenConf.<lambda> at 0x7f810efb7790>, as_digests=False)
```

List genomes and assets available remotely on all servers the object subscribes to
#### Parameters:

- `get_url` (`function(serverUrl, operationId) -> str`):  how to determineURL request, given server URL and endpoint operationID
- `genome` (`list[str] | str`):  genomes that the assets should be found for
- `order` (`function(str) -> object`):  how to key genome IDs and assetnames for sort


#### Returns:

- `dict[OrderedDict[list]]`:  remotely available genomes and assetskeyed by genome keyed by source server endpoint




```python
def plugins(self)
```

Plugins registered by entry points in the current Python env
#### Returns:

- `dict[dict[function(refgenconf.RefGenConf)]]`:  dict which keysare names of all possible hooks and values are dicts mapping registered functions names to their values




```python
def populate(self, glob)
```

Populates *local* refgenie references from refgenie://genome/asset.seek_key:tag registry paths
#### Parameters:

- `glob` (`dict | str | list`):  String which may contain refgenie registry paths asvalues; or a dict, for which values may contain refgenie registry paths. Dict include nested dicts.


#### Returns:

- `dict | str | list`:  modified input dict with refgenie paths populated




```python
def populater(self, glob, remote_class=None)
```

Populates *remote* refgenie references from refgenie://genome/asset:tag registry paths
#### Parameters:

- `glob` (`dict | str | list`):  String which may contain refgenie registry paths asvalues; or a dict, for which values may contain refgenie registry paths. Dict include nested dicts.
- `remote_class` (`str`):  remote data provider class, e.g. 'http' or 's3'


#### Returns:

- `dict | str | list`:  modified input dict with refgenie paths populated




```python
def pull(self, genome, asset, tag, unpack=True, force=None, force_large=None, size_cutoff=10, get_json_url=<function RefGenConf.<lambda> at 0x7f810efb7a60>, build_signal_handler=<function _handle_sigint at 0x7f810f4541f0>)
```

Download and possibly unpack one or more assets for a given ref gen.
#### Parameters:

- `genome` (`str`):  name of a reference genome assembly of interest
- `asset` (`str`):  name of particular asset to fetch
- `tag` (`str`):  name of particular tag to fetch
- `unpack` (`bool`):  whether to unpack a tarball
- `force` (`bool | NoneType`):  how to handle case in which asset pathalready exists; null for prompt (on a per-asset basis), False to effectively auto-reply No to the prompt to replace existing file, and True to auto-replay Yes for existing asset replacement.
- `force_large` (`bool | NoneType`):  how to handle case in large (> 5GB)asset is to be pulled; null for prompt (on a per-asset basis), False to effectively auto-reply No to the prompt, and True to auto-replay Yes
- `size_cutoff` (`float`):  maximum archive file size to download withno prompt
- `get_json_url` (`function(str, str) -> str`):  how to build URL fromgenome server URL base, genome, and asset
- `build_signal_handler` (`function(str) -> function`):  how to createa signal handler to use during the download; the single argument to this function factory is the download filepath


#### Returns:

- `(list[str], dict, str)`:  a list of genome, asset, tag namesand a key-value pair with which genome config file should be updated if pull succeeds, else asset key and a null value


#### Raises:

- `refgenconf.UnboundEnvironmentVariablesError`:  if genome folderpath contains any env. var. that's unbound
- `refgenconf.RefGenConfError`:  if the object update is requested ina non-writable state




```python
def remove(self, genome, asset, tag=None, relationships=True, files=True, force=False)
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
- `relationships` (`bool`):  whether the asset being removed shouldbe removed from its relatives as well
- `files` (`bool`):  whether the asset files from disk should be removed
- `force` (`bool`):  whether the removal prompts should be skipped


#### Returns:

- `RefGenConf`:  updated object


#### Raises:

- `TypeError`:  if genome argument type is not a list or str




```python
def remove_asset_from_relatives(self, genome, asset, tag)
```

Remove any relationship links associated with the selected asset
#### Parameters:

- `genome` (`str`):  genome to be removed from its relatives' relatives list
- `asset` (`str`):  asset to be removed from its relatives' relatives list
- `tag` (`str`):  tag to be removed from its relatives' relatives list




```python
def remove_genome_aliases(self, digest, aliases=None)
```

Remove alias for a specified genome digest. This method will remove the digest both from the genomes object and from the aliases mapping in tbe config
#### Parameters:

- `digest` (`str`):  genome digest to remove an alias for
- `aliases` (`list[str]`):  a collection to aliases to remove for thegenome. If not provided, all aliases for the digest will be remove


#### Returns:

- `bool`:  whether the removal has been performed




```python
def run_plugins(self, hook)
```

Runs all installed plugins for the specified hook.
#### Parameters:

- `hook` (`str`):  hook identifier




```python
def seek(self, genome_name, asset_name, tag_name=None, seek_key=None, strict_exists=None, enclosing_dir=False, all_aliases=False, check_exist=<function RefGenConf.<lambda> at 0x7f810efbcf70>)
```

Seek path to a specified genome-asset-tag alias
#### Parameters:

- `genome_name` (`str`):  name of a reference genome assembly of interest
- `asset_name` (`str`):  name of the particular asset to fetch
- `tag_name` (`str`):  name of the particular asset tag to fetch
- `seek_key` (`str`):  name of the particular subasset to fetch
- `strict_exists` (`bool | NoneType`):  how to handle case in whichpath doesn't exist; True to raise IOError, False to raise RuntimeWarning, and None to do nothing at all. Default: None (do not check).
- `check_exist` (`function(callable) -> bool`):  how to check forasset/path existence
- `enclosing_dir` (`bool`):  whether a path to the entire enclosingdirectory should be returned, e.g. for a fasta asset that has 3 seek_keys pointing to 3 files in an asset dir, that asset dir is returned
- `all_aliases` (`bool`):  whether to return paths to all asset aliases orjust the one for the specified 'genome_name` argument


#### Returns:

- `str`:  path to the asset


#### Raises:

- `TypeError`:  if the existence check is not a one-arg function
- `refgenconf.MissingGenomeError`:  if the named assembly isn't knownto this configuration instance
- `refgenconf.MissingAssetError`:  if the names assembly is known tothis configuration instance, but the requested asset is unknown




```python
def seek_src(self, genome_name, asset_name, tag_name=None, seek_key=None, strict_exists=None, enclosing_dir=False, check_exist=<function RefGenConf.<lambda> at 0x7f810efb71f0>)
```

Seek path to a specified genome-asset-tag
#### Parameters:

- `genome_name` (`str`):  name of a reference genome assembly of interest
- `asset_name` (`str`):  name of the particular asset to fetch
- `tag_name` (`str`):  name of the particular asset tag to fetch
- `seek_key` (`str`):  name of the particular subasset to fetch
- `strict_exists` (`bool | NoneType`):  how to handle case in whichpath doesn't exist; True to raise IOError, False to raise RuntimeWarning, and None to do nothing at all. Default: None (do not check).
- `check_exist` (`function(callable) -> bool`):  how to check forasset/path existence
- `enclosing_dir` (`bool`):  whether a path to the entire enclosingdirectory should be returned, e.g. for a fasta asset that has 3 seek_keys pointing to 3 files in an asset dir, that asset dir is returned


#### Returns:

- `str`:  path to the asset


#### Raises:

- `TypeError`:  if the existence check is not a one-arg function
- `refgenconf.MissingGenomeError`:  if the named assembly isn't knownto this configuration instance
- `refgenconf.MissingAssetError`:  if the names assembly is known tothis configuration instance, but the requested asset is unknown




```python
def seekr(self, genome_name, asset_name, tag_name=None, seek_key=None, remote_class='http', get_url=<function RefGenConf.<lambda> at 0x7f810efb70d0>)
```

Seek a remote path to a specified genome/asset.seek_key:tag
#### Parameters:

- `genome_name` (`str`):  name of a reference genome assembly of interest
- `asset_name` (`str`):  name of the particular asset to fetch
- `tag_name` (`str`):  name of the particular asset tag to fetch
- `seek_key` (`str`):  name of the particular subasset to fetch
- `remote_class` (`str`):  remote data provider class, e.g. 'http' or 's3'
- `get_url` (`function(serverUrl, operationId) -> str`):  how to determineURL request, given server URL and endpoint operationID


#### Returns:

- `str`:  path to the asset




```python
def set_default_pointer(self, genome, asset, tag, force_exists=False, force_digest=None, force_fasta=False)
```

Point to the selected tag by default
#### Parameters:

- `genome` (`str`):  name of a reference genome assembly of interest
- `asset` (`str`):  name of the particular asset of interest
- `tag` (`str`):  name of the particular asset tag to point to by default
- `force_digest` (`str`):  digest to force update of. The alias willnot be converted to the digest, even if provided.
- `force_fasta` (`bool`):  whether setting a default tag for a fasta assetshould be forced. Beware: This could lead to genome identity issues
- `force_exists` (`bool`):  whether the default tag change should beforced (even if it exists)




```python
def set_genome_alias(self, genome, digest=None, servers=None, overwrite=False, reset_digest=False, create_genome=False, no_write=False, get_json_url=<function RefGenConf.<lambda> at 0x7f810efb7d30>)
```

Assign a human-readable alias to a genome identifier.

Genomes are identified by a unique identifier which is derived from the
FASTA file (part of fasta asset). This way we can ensure genome
provenance and compatibility with the server. This function maps a
human-readable identifier to make referring to the genomes easier.
#### Parameters:

- `genome` (`str`):  name of the genome to assign to an identifier
- `digest` (`str`):  identifier to use
- `overwrite` (`bool`):  whether all the previously set aliases should beremoved and just the current one stored
- `no_write` (`bool`):  whether to skip writing the alias to the file


#### Returns:

- `bool`:  whether the alias has been established




```python
def subscribe(self, urls, reset=False, no_write=False)
```

Add URLs the list of genome_servers.

Use reset argument to overwrite the current list.
Otherwise the current one will be appended to.
#### Parameters:

- `urls` (`list[str] | str`):  urls to update the genome_servers list with
- `reset` (`bool`):  whether the current list should be overwritten




```python
def tag(self, genome, asset, tag, new_tag, files=True, force=False)
```

Retags the asset selected by the tag with the new_tag. Prompts if default already exists and overrides upon confirmation.

This method does not override the original asset entry in the RefGenConf
object. It creates its copy and tags it with the new_tag.
Additionally, if the retagged asset has any children their parent will
be retagged as new_tag that was introduced upon this method execution.
By default, the files on disk will be also renamed to reflect the
genome configuration file changes
#### Parameters:

- `genome` (`str`):  name of a reference genome assembly of interest
- `asset` (`str`):  name of particular asset of interest
- `tag` (`str`):  name of the tag that identifies the asset of interest
- `new_tag` (`str`):  name of particular the new tag
- `files` (`bool`):  whether the asset files on disk should be renamed


#### Returns:

- `bool`:  a logical indicating whether the tagging was successful


#### Raises:

- `ValueError`:  when the original tag is not specified




```python
def unsubscribe(self, urls, no_write=False)
```

Remove URLs the list of genome_servers.
#### Parameters:

- `urls` (`list[str] | str`):  urls to update the genome_servers list with




```python
def update_assets(self, genome, asset=None, data=None, force_digest=None)
```

Updates the genomes in RefGenConf object at any level. If a requested genome-asset mapping is missing, it will be created
#### Parameters:

- `genome` (`str`):  genome to be added/updated
- `asset` (`str`):  asset to be added/updated
- `force_digest` (`str`):  digest to force update of. The alias willnot be converted to the digest, even if provided.
- `data` (`Mapping`):  data to be added/updated


#### Returns:

- `RefGenConf`:  updated object




```python
def update_genomes(self, genome, data=None, force_digest=None)
```

Updates the genomes in RefGenConf object at any level. If a requested genome is missing, it will be added
#### Parameters:

- `genome` (`str`):  genome to be added/updated
- `force_digest` (`str`):  digest to force update of. The alias willnot be converted to the digest, even if provided.
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
- `data` (`list`):  asset parents or children to be added/updated
- `children` (`bool`):  a logical indicating whether the relationship to beadded is 'children'


#### Returns:

- `RefGenConf`:  updated object




```python
def update_seek_keys(self, genome, asset, tag=None, keys=None, force_digest=None)
```

A convenience method which wraps the updated assets and uses it to update the seek keys for a tagged asset.
#### Parameters:

- `genome` (`str`):  genome to be added/updated
- `asset` (`str`):  asset to be added/updated
- `tag` (`str`):  tag to be added/updated
- `force_digest` (`str`):  digest to force update of. The alias willnot be converted to the digest, even if provided.
- `keys` (`Mapping`):  seek_keys to be added/updated


#### Returns:

- `RefGenConf`:  updated object




```python
def update_tags(self, genome, asset=None, tag=None, data=None, force_digest=None)
```

Updates the genomes in RefGenConf object at any level. If a requested genome-asset-tag mapping is missing, it will be created
#### Parameters:

- `genome` (`str`):  genome to be added/updated
- `asset` (`str`):  asset to be added/updated
- `tag` (`str`):  tag to be added/updated
- `force_digest` (`str`):  digest to force update of. The alias willnot be converted to the digest, even if provided.
- `data` (`Mapping`):  data to be added/updated


#### Returns:

- `RefGenConf`:  updated object




```python
def writable(self)
```

Return writability flag or None if not set
#### Returns:

- `bool | None`:  whether the object is writable now




```python
def write(self, filepath=None)
```

Write the contents to a file. If pre- and post-update plugins are defined, they will be executed automatically
#### Parameters:

- `filepath` (`str`):  a file path to write to


#### Returns:

- `str`:  the path to the created files


#### Raises:

- `OSError`:  when the object has been created in a read only mode or otherprocess has locked the file
- `TypeError`:  when the filepath cannot be determined.This takes place only if YacAttMap initialized with a Mapping as an input, not read from file.
- `OSError`:  when the write is called on an object with no write capabilitiesor when writing to a file that is locked by a different object




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




```python
def get_dir_digest(path, pm=None)
```

Generate a MD5 digest that reflects just the contents of the files in the selected directory.
#### Parameters:

- `path` (`str`):  path to the directory to digest
- `pm` (`pypiper.PipelineManager`):  a pipeline object, optional.The subprocess module will be used if not provided


#### Returns:

- `str`:  a digest, e.g. a3c46f201a3ce7831d85cf4a125aa334




```python
def looper_refgenie_populate(namespaces)
```

A looper plugin that populates refgenie references in a PEP from refgenie://genome/asset:tag registry paths. This can be used to convert all refgenie references into their local paths at the looper stage, so the final paths are passed to the workflow. This way the workflow does not need to depend on refgenie to resolve the paths. This is useful for example for CWL pipelines, which are built to have paths resolved outside the workflow.

The namespaces structure required to run the plugin is:
`namespaces["pipeline"]["var_templates"]["refgenie_config"]`
#### Parameters:

- `namespaces` (`Mapping`):  a nested variable namespaces dict


#### Returns:

- `dict`:  sample namespace dict


#### Raises:

- `TypeError`:  if the input namespaces is not a mapping
- `KeyError`:  if the namespaces mapping does not include 'pipeline'
- `NotImplementedError`:  if 'var_templates' key is missing in the 'pipeline' namespace or'refgenie_config' is missing in 'var_templates' section.







*Version Information: `refgenconf` v0.12.2, generated by `lucidoc` v0.4.3*