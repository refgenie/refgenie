# Refgenie tutorial

The `Refgenie` class is a key component of the Refgenie package, which is used for managing and organizing reference genome files. It provides a set of methods and functionalities to interact with reference genome assets, and other related resources.

## Purpose of this file

This file serves as a tutorial for using the `Refgenie` Python API. It demonstrates how to set up a temporary directory for storing reference genome assets, configure the `Refgenie` instance, and perform various operations such as listing available assets, retrieving asset information, and managing data channels. 

In order to learn more about any of the contepts indicated in the code, please refer to a specific section of the documentation.

## Installation

Before the package gets released, clone the repository install, for example using `uv`:
  
```bash
git clone <repo_url>
cd refgenie1
uv pip install .
```

## Configuration

First, let's create a temporary directory that will be used to store the refgenie assets.



```python
from pathlib import Path
from rich import print
import os
```

Let's set a temporary directory to store the refgenie assets.


```python
from tempfile import TemporaryDirectory

REFGENIE_CODE_PATH = Path.cwd().parent / "refgenie"

# set the environment variable
archive_tmp_dir = TemporaryDirectory(prefix="refgenie_archive_demo_").name
os.environ["REFGENIE_GENOME_STAGE_FOLDER"] = archive_tmp_dir
tmp_dir = TemporaryDirectory(prefix="refgenie_demo_").name
os.environ["REFGENIE_GENOME_FOLDER"] = tmp_dir
# set the REFGENIE_DB_CONFIG_PATH to a sqlite config file in the refgenie package
os.environ["REFGENIE_DB_CONFIG_PATH"] = (REFGENIE_CODE_PATH / "config" / "sqlite_config.yaml").as_posix()

```

Let's inspect the refgenie configuration object.


```python
from refgenie.config import config

print(config)
```


<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #800080; text-decoration-color: #800080; font-weight: bold">RefgenieConfig</span><span style="font-weight: bold">(</span>
    <span style="color: #808000; text-decoration-color: #808000">log_level</span>=<span style="font-weight: bold">&lt;</span><span style="color: #ff00ff; text-decoration-color: #ff00ff; font-weight: bold">LogLevel.INFO:</span><span style="color: #000000; text-decoration-color: #000000"> </span><span style="color: #008000; text-decoration-color: #008000">'INFO'</span><span style="font-weight: bold">&gt;</span>,
    <span style="color: #808000; text-decoration-color: #808000">genome_folder</span>=<span style="color: #800080; text-decoration-color: #800080; font-weight: bold">PosixPath</span><span style="font-weight: bold">(</span><span style="color: #008000; text-decoration-color: #008000">'/tmp/refgenie_demo_2hzk45es'</span><span style="font-weight: bold">)</span>,
    <span style="color: #808000; text-decoration-color: #808000">genome_stage_folder</span>=<span style="color: #800080; text-decoration-color: #800080; font-weight: bold">PosixPath</span><span style="font-weight: bold">(</span><span style="color: #008000; text-decoration-color: #008000">'/tmp/refgenie_archive_demo_3ps6xfvv'</span><span style="font-weight: bold">)</span>,
    <span style="color: #808000; text-decoration-color: #808000">database_config_path</span>=<span style="color: #800080; text-decoration-color: #800080; font-weight: bold">PosixPath</span><span style="font-weight: bold">(</span><span style="color: #008000; text-decoration-color: #008000">'/home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/config/sqlite</span>
<span style="color: #008000; text-decoration-color: #008000">_config.yaml'</span><span style="font-weight: bold">)</span>,
    <span style="color: #808000; text-decoration-color: #808000">bridge_mode</span>=<span style="color: #008000; text-decoration-color: #008000">'read'</span>,
    <span style="color: #808000; text-decoration-color: #808000">bridge_origins</span>=<span style="color: #008000; text-decoration-color: #008000">'https://refgenie.org,https://ui.refgenie.org'</span>,
    <span style="color: #808000; text-decoration-color: #808000">bridge_origin_regex</span>=<span style="color: #008000; text-decoration-color: #008000">''</span>,
    <span style="color: #808000; text-decoration-color: #808000">bridge_expose_paths</span>=<span style="color: #ff0000; text-decoration-color: #ff0000; font-style: italic">False</span>
<span style="font-weight: bold">)</span>
</pre>



### Database backend

As you can see, refgenie configuration points to a database configuration file, as by default refgenie is backed by a SQLite database.

Let's inspect the refgenie database configuration file.


```python
%cat {config.database_config_path}
```

    type: sqlite
    path: ~/refgenie_db/refgenie


Make sure the directory where the SQLite database file is stored exists, and create it if it doesn't.


```python
!refgenie purge --force
!rm -rf ~/refgenie_db
!mkdir -p ~/refgenie_db
```

    Traceback (most recent call last):
      File "/home/nsheff/.local/bin/refgenie", line 5, in <module>
        from refgenie.cli.cli_pydantic import main_cli
    ModuleNotFoundError: No module named 'refgenie.cli.cli_pydantic'


In practice, you don't even need to create the configuration file manually, as refgenie ships with a default configuration file that is used if no configuration file is provided. Just as we've seen above.

For production deployments you may want to use a different database backend, such as MySQL or PostgreSQL. In this case, you can provide the database configuration file path by setting `REFGENIE_DB_CONFIG_PATH` environment variable, or even set/override the database engine using `database_engine` in the `Refgenie` constructor. The object must be a `sqlalchemy.engine.Engine` object.



### Refgenieserver client

Similarly, refgenie ships with a Refgenieserver client, which is used by default to retrieve remote genome assets and does not need to be replaced in majority of use cases. However, you can provide a custom URL-client mapping to `Refgenie` constructor, by setting the `server_client_mapping` argument. Please note that, the clients need to follow a specific interface, defined in `refgenie.server.ServerClient` protocol. More details below.


```python
from refgenie.managers.sources.client import ServerClient
from rich import inspect

inspect(
    ServerClient,
    methods=True,
    docs=True,
    help=True,
    title="ServerClient Protocol structure",
)
```


<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">╭──────────────────────────────────────── ServerClient Protocol structure ────────────────────────────────────────╮</span>
<span style="color: #000080; text-decoration-color: #000080">│</span> <span style="color: #00ffff; text-decoration-color: #00ffff; font-style: italic">class </span><span style="color: #800000; text-decoration-color: #800000; font-weight: bold">ServerClient</span><span style="font-weight: bold">(</span>*args, **kwargs<span style="font-weight: bold">)</span>:                                                                            <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>                                                                                                                 <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span> <span style="color: #008080; text-decoration-color: #008080">Protocol for the server client.</span>                                                                                 <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>                                                                                                                 <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>             <span style="color: #808000; text-decoration-color: #808000; font-style: italic">server_url</span> = <span style="font-weight: bold">&lt;</span><span style="color: #ff00ff; text-decoration-color: #ff00ff; font-weight: bold">property</span><span style="color: #000000; text-decoration-color: #000000"> object at </span><span style="color: #008080; text-decoration-color: #008080; font-weight: bold">0x7ad9be7b03b0</span><span style="font-weight: bold">&gt;</span>                                                    <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>          <span style="color: #808000; text-decoration-color: #808000; font-style: italic">download_file</span> = <span style="color: #00ffff; text-decoration-color: #00ffff; font-style: italic">def </span><span style="color: #800000; text-decoration-color: #800000; font-weight: bold">download_file</span><span style="font-weight: bold">(</span>self, asset_digest: str, file_path: str, output_path: pathlib.Path<span style="font-weight: bold">)</span>  <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>                          -&gt; pathlib.Path: <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">Download a single file from a file-mode asset.</span>                        <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span> <span style="color: #808000; text-decoration-color: #808000; font-style: italic">download_with_progress</span> = <span style="color: #00ffff; text-decoration-color: #00ffff; font-style: italic">def </span><span style="color: #800000; text-decoration-color: #800000; font-weight: bold">download_with_progress</span><span style="font-weight: bold">(</span>self, operation_id: str, output_path: pathlib.Path, params: <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>                          dict | <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span> = <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span>, url_format_params: dict<span style="font-weight: bold">[</span>str, str | <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span><span style="font-weight: bold">]</span> = <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span>, name: str | <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span>  <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>                          = <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span><span style="font-weight: bold">)</span> -&gt; pathlib.Path:                                                               <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">Download the asset served by the given operation_id to output_path,</span>                    <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">showing progress along the way.</span>                                                        <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>                                                                                                                 <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">Implementations should report byte progress to `refgenie.progress`</span>                     <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">whenever `</span><span style="color: #bf7fbf; text-decoration-color: #bf7fbf; font-weight: bold">refgenie.progress.active_sink</span><span style="color: #7f7f7f; text-decoration-color: #7f7f7f; font-weight: bold">()</span><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">` is not </span><span style="color: #bf7fbf; text-decoration-color: #bf7fbf; font-style: italic">None</span><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">, and should not</span>                 <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">build a `rich` live display in that case -- something other than a</span>                     <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">terminal is consuming the progress.</span>                                                    <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>                    <span style="color: #808000; text-decoration-color: #808000; font-style: italic">get</span> = <span style="color: #00ffff; text-decoration-color: #00ffff; font-style: italic">def </span><span style="color: #800000; text-decoration-color: #800000; font-weight: bold">get</span><span style="font-weight: bold">(</span>self, operation_id: str, params: dict | <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span> = <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span>, url_format_params:        <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>                          dict<span style="font-weight: bold">[</span>str, str | <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span><span style="font-weight: bold">]</span> = <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span><span style="font-weight: bold">)</span> -&gt; dict: <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">Send a GET request to the specified operation </span>  <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">ID.</span>                                                                                    <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>        <span style="color: #808000; text-decoration-color: #808000; font-style: italic">get_all_aliases</span> = <span style="color: #00ffff; text-decoration-color: #00ffff; font-style: italic">def </span><span style="color: #800000; text-decoration-color: #800000; font-weight: bold">get_all_aliases</span><span style="font-weight: bold">(</span>self, params: dict | <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span> = <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span><span style="font-weight: bold">)</span> -&gt; list<span style="font-weight: bold">[</span>dict<span style="font-weight: bold">]</span>: <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">Get all aliases </span>  <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">using pagination.</span>                                                                      <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>   <span style="color: #808000; text-decoration-color: #808000; font-style: italic">get_all_asset_groups</span> = <span style="color: #00ffff; text-decoration-color: #00ffff; font-style: italic">def </span><span style="color: #800000; text-decoration-color: #800000; font-weight: bold">get_all_asset_groups</span><span style="font-weight: bold">(</span>self, params: dict | <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span> = <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span><span style="font-weight: bold">)</span> -&gt; list<span style="font-weight: bold">[</span>dict<span style="font-weight: bold">]</span>: <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">Get all </span>     <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset groups using pagination.</span>                                                         <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>         <span style="color: #808000; text-decoration-color: #808000; font-style: italic">get_all_assets</span> = <span style="color: #00ffff; text-decoration-color: #00ffff; font-style: italic">def </span><span style="color: #800000; text-decoration-color: #800000; font-weight: bold">get_all_assets</span><span style="font-weight: bold">(</span>self, params: dict | <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span> = <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span><span style="font-weight: bold">)</span> -&gt; list<span style="font-weight: bold">[</span>dict<span style="font-weight: bold">]</span>: <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">Get all assets </span>    <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">using pagination.</span>                                                                      <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>        <span style="color: #808000; text-decoration-color: #808000; font-style: italic">get_all_genomes</span> = <span style="color: #00ffff; text-decoration-color: #00ffff; font-style: italic">def </span><span style="color: #800000; text-decoration-color: #800000; font-weight: bold">get_all_genomes</span><span style="font-weight: bold">(</span>self, params: dict | <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span> = <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span><span style="font-weight: bold">)</span> -&gt; list<span style="font-weight: bold">[</span>dict<span style="font-weight: bold">]</span>: <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">Get all genomes </span>  <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">using pagination.</span>                                                                      <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>  <span style="color: #808000; text-decoration-color: #808000; font-style: italic">get_all_staged_assets</span> = <span style="color: #00ffff; text-decoration-color: #00ffff; font-style: italic">def </span><span style="color: #800000; text-decoration-color: #800000; font-weight: bold">get_all_staged_assets</span><span style="font-weight: bold">(</span>self, params: dict | <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span> = <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span><span style="font-weight: bold">)</span> -&gt; list<span style="font-weight: bold">[</span>dict<span style="font-weight: bold">]</span>: <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">Get all </span>    <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">staged assets using pagination.</span>                                                        <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>    <span style="color: #808000; text-decoration-color: #808000; font-style: italic">get_asset_file_list</span> = <span style="color: #00ffff; text-decoration-color: #00ffff; font-style: italic">def </span><span style="color: #800000; text-decoration-color: #800000; font-weight: bold">get_asset_file_list</span><span style="font-weight: bold">(</span>self, asset_digest: str<span style="font-weight: bold">)</span> -&gt; list<span style="font-weight: bold">[</span>str<span style="font-weight: bold">]</span>: <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">Get the list of files </span>  <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">available for a file-mode asset.</span>                                                       <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>       <span style="color: #808000; text-decoration-color: #808000; font-style: italic">get_asset_groups</span> = <span style="color: #00ffff; text-decoration-color: #00ffff; font-style: italic">def </span><span style="color: #800000; text-decoration-color: #800000; font-weight: bold">get_asset_groups</span><span style="font-weight: bold">(</span>self, params: dict | <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span> = <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span><span style="font-weight: bold">)</span> -&gt; list<span style="font-weight: bold">[</span>dict<span style="font-weight: bold">]</span>: <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">Get one page of </span> <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset groups, unwrapping the paginated response envelope.</span>                              <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>             <span style="color: #808000; text-decoration-color: #808000; font-style: italic">get_assets</span> = <span style="color: #00ffff; text-decoration-color: #00ffff; font-style: italic">def </span><span style="color: #800000; text-decoration-color: #800000; font-weight: bold">get_assets</span><span style="font-weight: bold">(</span>self, params: dict | <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span> = <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span><span style="font-weight: bold">)</span> -&gt; list<span style="font-weight: bold">[</span>dict<span style="font-weight: bold">]</span>: <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">Get one page of </span>       <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">assets, unwrapping the paginated response envelope.</span>                                    <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>      <span style="color: #808000; text-decoration-color: #808000; font-style: italic">get_staged_assets</span> = <span style="color: #00ffff; text-decoration-color: #00ffff; font-style: italic">def </span><span style="color: #800000; text-decoration-color: #800000; font-weight: bold">get_staged_assets</span><span style="font-weight: bold">(</span>self, params: dict | <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span> = <span style="color: #800080; text-decoration-color: #800080; font-style: italic">None</span><span style="font-weight: bold">)</span> -&gt; list<span style="font-weight: bold">[</span>dict<span style="font-weight: bold">]</span>: <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">Get one page of</span> <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">│</span>                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">staged assets, unwrapping the paginated response envelope.</span>                             <span style="color: #000080; text-decoration-color: #000080">│</span>
<span style="color: #000080; text-decoration-color: #000080">╰─────────────────────────────────────────────────────────────────────────────────────────────────────────────────╯</span>
</pre>



First, let's import the `Refgenie` class from the `refgenie` package.



```python
from refgenie import GenomeAlias, Refgenie

refgenie = Refgenie(suppress_migrations=True)
```

Let's ensure we start with a clean slate by removing any existing refgenie metadata and initializing a new refgenie instance.



```python
refgenie.init()  # initialize new refgenie instance
```


<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Genome folder ready: <span style="color: #800080; text-decoration-color: #800080">/tmp/</span><span style="color: #ff00ff; text-decoration-color: #ff00ff">refgenie_demo_2hzk45es</span>                                          <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/core/lifecycle.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">lifecycle.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/core/lifecycle.py#182" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">182</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Genome stage folder ready: <span style="color: #800080; text-decoration-color: #800080">/tmp/</span><span style="color: #ff00ff; text-decoration-color: #ff00ff">refgenie_archive_demo_3ps6xfvv</span>                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/core/lifecycle.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">lifecycle.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/core/lifecycle.py#186" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">186</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Initialized refgenie backend: <span style="color: #008000; text-decoration-color: #008000">'sqlite:////home/nsheff/refgenie_db/refgenie'</span>               <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/core/lifecycle.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">lifecycle.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/core/lifecycle.py#106" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">106</span></a>
</pre>



Before we can pull or build a `fasta` asset, refgenie needs to know the `fasta` asset class and recipe. Let's register them from the recipes repository (we'll add a whole data channel of recipes later, but for now we just need `fasta`):


```python
refgenie.asset_class.add(
    "https://github.com/refgenie/recipes/raw/refgenie1/asset_classes/fasta_asset_class.yaml"
)
refgenie.recipe.add(
    "https://github.com/refgenie/recipes/raw/refgenie1/recipes/fasta_asset_recipe.yaml"
)
```


<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://github.com/refgenie/recipes/raw/refgenie1/asset_classes/fasta_asset_class.yaml</span>            <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://github.com/refgenie/recipes/raw/refgenie1/asset_classes/fasta_asset_class.yaml</span>     <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 302 Found"</span>                                                                       <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://raw.githubusercontent.com/refgenie/recipes/refgenie1/asset_classes/fasta_asset_cla</span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">ss.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                                  <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'fasta'</span> asset class                                                           <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://github.com/refgenie/recipes/raw/refgenie1/recipes/fasta_asset_recipe.yaml</span>                 <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://github.com/refgenie/recipes/raw/refgenie1/recipes/fasta_asset_recipe.yaml</span>          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 302 Found"</span>                                                                       <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://raw.githubusercontent.com/refgenie/recipes/refgenie1/recipes/fasta_asset_recipe.ya</span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">ml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                                       <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'fasta'</span> recipe                                                                    <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>






    Recipe(id=1, name='fasta', version='0.1.0', description='DNA sequences in the FASTA format, exported from RefgetStore. Includes FASTA index (.fai) and chromosome sizes file. Requires genome to be initialized in RefgetStore first (via refgenie genome init).', output_asset_class_id=1, command_templates=['refgenie-build-fasta {{values.refget_store_path}} {{values.genome_digest}} {{values.output_folder}}'], input_params=None, input_files=None, input_assets=None, docker_image=None, custom_seek_keys={}, default_asset='default', inherent=None, updated_at=datetime.datetime(2026, 8, 21, 1, 55, 43, 254318, tzinfo=datetime.timezone.utc), created_at=datetime.datetime(2026, 8, 21, 1, 55, 43, 254321, tzinfo=datetime.timezone.utc))



Let's subscribe to the default refgenie server. This method will reach out to the server at the provided URL and query the OpenAPI specification to determine whether ther server is refgenie-compatible. If it is, the server will be added to the list of subscribed servers.

Note: there's currently no public compatible refgenieserver instance deployed, so the following code snippets use a local refgenieserver instance serving the latest API.



```python
refgenie.servers.subscribe("http://localhost:8000")
```


<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Subscribed to servers: <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">http://localhost:8000</span>                                           <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/configuration.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">configuration.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/configuration.py#55" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">55</span></a>
</pre>



And that's it! We have now configured a refgenie instance and subscribed to a refgenie-compatible server. We can now start using the refgenie instance to manage reference genome assets.

### Pull an asset

Let's initialize a new genome by pulling an asset of fasta class. This will create a new directory in the `data` subdirectory of the `genome_folder` and mirror it in the `alias` directory with symbolic links, rather than copies of the files.



```python
refgenie.pull(genome=GenomeAlias("rCRSd"), asset_group_name="fasta")
```


<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span> <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">http://localhost:8000/openapi.json</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                     <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Connected to server: <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">http://localhost:8000</span> <span style="color: #808000; text-decoration-color: #808000">title</span>=<span style="color: #008000; text-decoration-color: #008000">'Refgenieserver REST API'</span> <span style="color: #808000; text-decoration-color: #808000">version</span>=<span style="color: #008000; text-decoration-color: #008000">'1.0.0a1'</span> <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/sources/client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/sources/client.py#221" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">221</span></a>
         <span style="color: #808000; text-decoration-color: #808000">description</span>=<span style="color: #008000; text-decoration-color: #008000">'a web interface and RESTful API for reference genome assets'</span>                    <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">             </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #808000; text-decoration-color: #808000">WARNING </span> No local digest for genome alias: rCRSd. Setting genome identity with server:                <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/puller.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">puller.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/puller.py#452" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">452</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">http://localhost:8000</span>                                                                        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">             </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span> <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">http://localhost:8000/rgstore.json</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                     <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">http://localhost:8000</span> is not a RefgetStore root; looking for a server service-info.         <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/sources/genomes.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">genomes.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/sources/genomes.py#355" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">355</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span> <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">http://localhost:8000/service-info</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                     <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span> <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">http://localhost:8000/v4/aliases/rCRSd</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                 <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span> <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">http://localhost:8000/v4/genomes/jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP</span>        <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Resolved alias <span style="color: #008000; text-decoration-color: #008000">'rCRSd'</span> via server client: jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP         <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/genome_bootstrap.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">genome_bootstrap.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/genome_bootstrap.py#110" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">110</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Determined digest for rCRSd: jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP                      <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/genome_bootstrap.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">genome_bootstrap.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/genome_bootstrap.py#118" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">118</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span> <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">http://localhost:8000/rgstore.json</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                     <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">http://localhost:8000</span> is not a RefgetStore root; looking for a server service-info.         <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/sources/genomes.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">genomes.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/sources/genomes.py#355" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">355</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span> <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">http://localhost:8000/service-info</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                     <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Added genome: jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP                                                <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/genome.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">genome.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/genome.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Added alias: rCRSd                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/alias.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">alias.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/alias.py#305" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">305</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Querying server <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">http://localhost:8000</span> for jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP/fasta             <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/puller.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">puller.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/puller.py#691" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">691</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">http://localhost:8000/v4/asset_groups?genome_digest=jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP&amp;asset</span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">_group_name=fasta</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span> <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">http://localhost:8000/v4/assets?asset_group_id=1&amp;name=</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span> <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">http://localhost:8000/v4/relationships/03aa63e8accae38cf4d6bcbb7865ef85f95a54020c4f4d48921</span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">e67ce382f188c?expand=true</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">http://localhost:8000/v4/staged_assets?asset_digest=03aa63e8accae38cf4d6bcbb7865ef85f95a54</span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">020c4f4d48921e67ce382f188c</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                               <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace">/home/nsheff/.local/lib/python3.12/site-packages/rich/live.py:260: UserWarning: install "ipywidgets" for Jupyter 
support
  warnings.warn('install "ipywidgets" for Jupyter support')
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">http://localhost:8000/v4/archives/03aa63e8accae38cf4d6bcbb7865ef85f95a54020c4f4d48921e67ce</span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">382f188c/download</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"></pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Extracting asset tarball:                                                                   <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/puller.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">puller.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/puller.py#1014" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1014</span></a>
         <span style="color: #800080; text-decoration-color: #800080">/tmp/refgenie_demo_2hzk45es/data/jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP/fasta/</span><span style="color: #ff00ff; text-decoration-color: #ff00ff">fasta__default.tgz</span>  <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">              </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span> <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">http://localhost:8000/v4/genomes/jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP</span>        <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Initializing genome from FASTA file:                                                         <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/genome.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">genome.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/genome.py#417" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">417</span></a>
         <span style="color: #800080; text-decoration-color: #800080">/tmp/refgenie_demo_2hzk45es/data/jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP/fasta/03aa63e8accae38cf4d6</span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">             </span>
         <span style="color: #800080; text-decoration-color: #800080">bcbb7865ef85f95a54020c4f4d48921e67ce382f188c/</span><span style="color: #ff00ff; text-decoration-color: #ff00ff">jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP.fa</span>             <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">             </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Added: <span style="color: #008000; text-decoration-color: #008000">'rCRSd/fasta:default'</span>                                                                <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/content.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">content.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/content.py#381" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">381</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Set default asset: <span style="color: #008000; text-decoration-color: #008000">'jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP/fasta:default'</span>                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/seek.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">seek.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/seek.py#612" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">612</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Created alias directories: <span style="color: #800080; text-decoration-color: #800080">/tmp/refgenie_demo_2hzk45es/alias/rCRSd/fasta/</span><span style="color: #ff00ff; text-decoration-color: #ff00ff">default</span>           <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/symlinks.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">symlinks.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/symlinks.py#223" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">223</span></a>
</pre>






    Asset(name='default', description=None, size=34113, colocate=None, updated_at=datetime.datetime(2026, 8, 21, 1, 55, 43, 893614), path='data/jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP/fasta/03aa63e8accae38cf4d6bcbb7865ef85f95a54020c4f4d48921e67ce382f188c', digest='03aa63e8accae38cf4d6bcbb7865ef85f95a54020c4f4d48921e67ce382f188c', recipe_id=None, serving_modes_override=None, asset_group_id=1, created_at=datetime.datetime(2026, 8, 21, 1, 55, 43, 893622))



As you can see above, the genome has been initialized and `fasta` asset was pulled. Let's inspect the initialized genome.



```python
print(refgenie.genome.table())
```


<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="font-style: italic">                                       Genomes                                       </span>
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━┓
┃<span style="font-weight: bold"> Digest                           </span>┃<span style="font-weight: bold"> Aliases </span>┃<span style="font-weight: bold"> Source </span>┃<span style="font-weight: bold"> Species </span>┃<span style="font-weight: bold"> Description       </span>┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━┩
│ jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP │ rCRSd   │ local  │         │ rCRSd test genome │
└──────────────────────────────────┴─────────┴────────┴─────────┴───────────────────┘
</pre>



Now, that a `fasta` asset has been built for the `dm6` genome, let's add some custom asset classes and recipes to build an asset based on that.

### Add `bowtie2_index` asset class and recipe

By supplying a URL (`str` object) rather than a local path (`pathlib.Path` object), refgenie will grab the remote file and register it as if it was a local file.



```python
refgenie.asset_class.add(
    "https://github.com/refgenie/recipes/raw/refgenie1/asset_classes/bowtie2_index_asset_class.yaml"
)
refgenie.recipe.add(
    "https://github.com/refgenie/recipes/raw/refgenie1/recipes/bowtie2_index_asset_recipe.yaml"
)
```


<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://github.com/refgenie/recipes/raw/refgenie1/asset_classes/bowtie2_index_asset_class.yaml</span>    <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://github.com/refgenie/recipes/raw/refgenie1/asset_classes/bowtie2_index_asset_class.</span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 302 Found"</span>                                                                  <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://raw.githubusercontent.com/refgenie/recipes/refgenie1/asset_classes/bowtie2_index_a</span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">sset_class.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'bowtie2_index'</span> asset class                                                   <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://github.com/refgenie/recipes/raw/refgenie1/recipes/bowtie2_index_asset_recipe.yaml</span>         <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://github.com/refgenie/recipes/raw/refgenie1/recipes/bowtie2_index_asset_recipe.yaml</span>  <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 302 Found"</span>                                                                       <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://raw.githubusercontent.com/refgenie/recipes/refgenie1/recipes/bowtie2_index_asset_r</span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">ecipe.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                               <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'bowtie2_index'</span> recipe                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>






    Recipe(id=2, name='bowtie2_index', version='0.0.1', description='Genome index for bowtie2, produced with bowtie2-build', output_asset_class_id=2, command_templates=['bowtie2-build --threads {{values.params["threads"]}} {{values.genome_folder}}/{{values.assets["fasta"].seek_keys_dict["fasta"]}} {{values.output_folder}}/{{values.genome_digest}}'], input_params={'threads': {'description': 'Number of threads to use', 'default': 1}}, input_files=None, input_assets={'fasta': {'asset_class': 'fasta', 'description': 'fasta asset for genome', 'default': 'fasta'}}, docker_image='docker.io/databio/refgenie', custom_seek_keys={'version': "bowtie2-build --version | awk 'NR==1{print $3}'"}, default_asset='{{values.custom_seek_keys.version}}', inherent=None, updated_at=datetime.datetime(2026, 8, 21, 1, 55, 44, 426089, tzinfo=datetime.timezone.utc), created_at=datetime.datetime(2026, 8, 21, 1, 55, 44, 426090, tzinfo=datetime.timezone.utc))



Let's verify that it worked by listing the available asset classes and recipes:



```python
from rich import print

print(refgenie.recipe.table())
print(refgenie.asset_class.table())
```


<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="font-style: italic">                                                      Recipes                                                      </span>
┏━━━━━━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━┓
┃<span style="font-weight: bold">               </span>┃<span style="font-weight: bold">         </span>┃<span style="font-weight: bold"> Output asset     </span>┃<span style="font-weight: bold"> Input asset     </span>┃<span style="font-weight: bold">             </span>┃<span style="font-weight: bold">                  </span>┃<span style="font-weight: bold">                 </span>┃
┃<span style="font-weight: bold"> Name          </span>┃<span style="font-weight: bold"> Version </span>┃<span style="font-weight: bold"> class            </span>┃<span style="font-weight: bold"> classes         </span>┃<span style="font-weight: bold"> Input files </span>┃<span style="font-weight: bold"> Input params     </span>┃<span style="font-weight: bold"> Docker image    </span>┃
┡━━━━━━━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━┩
│ fasta         │ 0.1.0   │ fasta            │ <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">None</span>            │ <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">None</span>        │ <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">None</span>             │                 │
├───────────────┼─────────┼──────────────────┼─────────────────┼─────────────┼──────────────────┼─────────────────┤
│ bowtie2_index │ 0.0.1   │ bowtie2_index    │ • fasta <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">(fasta </span> │ <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">None</span>        │ • threads        │ docker.io/data… │
│               │         │                  │ <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset for </span>      │             │ <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">(Number of </span>      │                 │
│               │         │                  │ <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">genome) </span>        │             │ <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">threads to use) </span> │                 │
│               │         │                  │ <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">default=fasta</span>   │             │ <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">default=1</span>        │                 │
└───────────────┴─────────┴──────────────────┴─────────────────┴─────────────┴──────────────────┴─────────────────┘
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="font-style: italic">                                                   Asset Classes                                                   </span>
┏━━━━━━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃<span style="font-weight: bold"> Name          </span>┃<span style="font-weight: bold"> Version </span>┃<span style="font-weight: bold"> Seek keys               </span>┃<span style="font-weight: bold"> Description                                                 </span>┃
┡━━━━━━━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ fasta         │ 0.1.0   │ fasta, fai, chrom_sizes │ Sequences in the FASTA format, indexed FASTA and chromosome │
│               │         │                         │ sizes file                                                  │
├───────────────┼─────────┼─────────────────────────┼─────────────────────────────────────────────────────────────┤
│ bowtie2_index │ 0.0.1   │ bowtie2_index           │ Genome index for bowtie2, produced with bowtie2-build       │
└───────────────┴─────────┴─────────────────────────┴─────────────────────────────────────────────────────────────┘
</pre>



## Build a `fasta` asset



```python
from refgenie import BuildParams
from pathlib import Path

# Step 1: Initialize genome from FASTA
refgenie.genome.initialize_genome(
    fasta_file_path=REFGENIE_CODE_PATH.parent / "tests/data/t7.fa",
    alias_names=["t7"],
    description="Genome of T7 phage",
)

# Step 2: Build fasta asset (exports from RefgetStore, no input file needed)
refgenie.build_asset(
    recipe_name="fasta",
    genome_name="t7",
    asset_group_name="fasta",
)
```


<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Initializing genome from FASTA file:                                                         <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/genome.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">genome.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/genome.py#417" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">417</span></a>
         <span style="color: #800080; text-decoration-color: #800080">/home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/tests/data/</span><span style="color: #ff00ff; text-decoration-color: #ff00ff">t7.fa</span>                    <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">             </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Added genome: kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB                                                <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/genome.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">genome.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/genome.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Added alias: t7                                                                               <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/alias.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">alias.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/alias.py#305" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">305</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Building <span style="color: #008000; text-decoration-color: #008000">'t7/fasta'</span> using recipe <span style="color: #008000; text-decoration-color: #008000">'fasta (v0.1.0)'</span>                                           <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/builder.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">builder.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/builder.py#513" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">513</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Could not locate config file.                                                                <a href="file:///home/nsheff/Dropbox/workspaces/pepkit/repos/yacman/yacman/yacman.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">yacman.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/pepkit/repos/yacman/yacman/yacman.py#674" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">674</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> No output schema supplied. Running in schema-optional mode <span style="font-weight: bold">(</span><span style="color: #808000; text-decoration-color: #808000">validate_results</span>=<span style="color: #ff0000; text-decoration-color: #ff0000; font-style: italic">False</span><span style="font-weight: bold">)</span>.      <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/pipestat/pipestat.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">pipestat.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/pipestat/pipestat.py#1008" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1008</span></a>
         Results will not be validated against a schema.                                           <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">                </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Initializing results file                                                               <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/pipestat/backends/file_backend/filebackend.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">filebackend.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/pipestat/backends/file_backend/filebackend.py#785" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">785</span></a>
         <span style="color: #008000; text-decoration-color: #008000">'/tmp/refgenie_demo_2hzk45es/builds/t7/fasta/default/stats.yaml'</span>                        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">                  </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> File does not exist, but create_file is true. Creating<span style="color: #808000; text-decoration-color: #808000">...</span>                              <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/ubiquerg/file_locking.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">file_locking.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/ubiquerg/file_locking.py#313" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">313</span></a>
</pre>



    ### Pipeline run code and environment:
    


    *          Command: `/home/nsheff/.local/lib/python3.12/site-packages/ipykernel_launcher.py -f /tmp/tmp30pcozza.json --HistoryManager.hist_file=:memory:`


    *     Compute host: `zither`


    *      Working dir: `/home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/docs`


    *        Outfolder: `/tmp/refgenie_demo_2hzk45es/builds/t7/fasta/default/`


    *         Log file: `/tmp/refgenie_demo_2hzk45es/builds/t7/fasta/default/refgenie_t7_fasta_default_log.md`


    *       Start time:  (08-20 21:55:44) elapsed: 0.0 _TIME_


    
    ### Version log:
    


    *   Python version: `3.12.3`


    *      Pypiper dir: `/home/nsheff/.local/lib/python3.12/site-packages/pypiper`


    *  Pypiper version: `0.15.1`


    * Pipestat version: `0.13.1`


    *     Pipeline dir: `/home/nsheff/.local/lib/python3.12/site-packages`


    * Pipeline version: 


    
    ### Arguments passed to pipeline:
    


    
    ### Initialized Pipestat Object:
    


    * PipestatManager (refgenie_t7_fasta_default)


    * Backend: File


    *  - results: /tmp/refgenie_demo_2hzk45es/builds/t7/fasta/default/stats.yaml


    *  - status: /tmp/refgenie_demo_2hzk45es/builds/t7/fasta/default


    * Multiple Pipelines Allowed: False


    * Pipeline name: refgenie_t7_fasta_default


    * Pipeline type: sample


    * Status Schema key: None


    * Results formatter: default_formatter


    * Results schema source: None


    * Status schema source: None


    * Records count: 0


    * Sample name: DEFAULT_SAMPLE_NAME
    


    
    ----------------------------------------
    


    Target to produce: `/tmp/refgenie_demo_2hzk45es/builds/t7/fasta/default/t7_fasta__default.flag`  


    
    > `refgenie-build-fasta /tmp/refgenie_demo_2hzk45es/.refget_store kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB /tmp/refgenie_demo_2hzk45es/data/kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB/fasta/default` (1045421)


    <pre>


    Loading collection metadata kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB...


    Loading sequence k0lfTnVEgjj6zguTU2MWdfYLenQFT6zW...


    </pre>


    Command completed. Elapsed time: 0:00:00. Running peak memory: 0.015GB.  
      PID: 1045421;	Command: refgenie-build-fasta;	Return code: 0;	Memory used: 0.015GB
    



<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Asset <span style="color: #008000; text-decoration-color: #008000">'t7/fasta:default'</span> build succeeded                                                    <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/builder.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">builder.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/builder.py#666" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">666</span></a>
</pre>



    
    ### Pipeline completed. Epilogue


    *        Elapsed time (this run):  0:00:00


    *  Total elapsed time (all runs):  0:00:00


    *         Peak memory (this run):  0.0154 GB


    *        Pipeline completed time: 2026-08-20 21:55:44



<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"></pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Added: <span style="color: #008000; text-decoration-color: #008000">'t7/fasta:default'</span>                                                                   <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/content.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">content.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/content.py#381" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">381</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Set default asset: <span style="color: #008000; text-decoration-color: #008000">'kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB/fasta:default'</span>                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/seek.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">seek.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/seek.py#612" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">612</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Added asset: <span style="color: #008000; text-decoration-color: #008000">'kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB/fasta:default'</span>                               <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/builder.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">builder.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/builder.py#705" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">705</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Created alias directories: <span style="color: #800080; text-decoration-color: #800080">/tmp/refgenie_demo_2hzk45es/alias/t7/fasta/</span><span style="color: #ff00ff; text-decoration-color: #ff00ff">default</span>              <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/symlinks.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">symlinks.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/symlinks.py#223" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">223</span></a>
</pre>






    Asset(name='default', description='DNA sequences in the FASTA format, exported from RefgetStore. Includes FASTA index (.fai) and chromosome sizes file. Requires genome to be initialized in RefgetStore first (via refgenie genome init).', size=40486, colocate=None, updated_at=datetime.datetime(2026, 8, 21, 1, 55, 44, 814061), path='data/kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB/fasta/949c705a1c1668f0ea3d7da883d1a65e902b945b67de576407349ff72a8b1bf1', digest='949c705a1c1668f0ea3d7da883d1a65e902b945b67de576407349ff72a8b1bf1', recipe_id=1, serving_modes_override=None, asset_group_id=2, created_at=datetime.datetime(2026, 8, 21, 1, 55, 44, 814074))



### Build a `bowtie2_index` asset

The `bowtie2_index` asset class and recipe have been added successfully. Let's build the `bowtie2_index` asset for the `dm6` genome.



```python
from refgenie.models import BuildParams

refgenie.build_asset(
    recipe_name="bowtie2_index",
    genome_name="t7",
    asset_group_name="bowtie2_index",
    params=BuildParams(params={"threads": 8}),
    stage=True,  # stage the asset right after building
)
```


<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Building <span style="color: #008000; text-decoration-color: #008000">'t7/bowtie2_index'</span> using recipe <span style="color: #008000; text-decoration-color: #008000">'bowtie2_index (v0.0.1)'</span>                           <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/builder.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">builder.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/builder.py#513" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">513</span></a>
</pre>



    [90m[[0m2026-08-21T01:55:44Z [33mWARN [0m bulker::shimlink[90m][0m 'bowtie2-build': `docker_command` is deprecated; use the `entrypoint` manifest field instead



<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Could not locate config file.                                                                <a href="file:///home/nsheff/Dropbox/workspaces/pepkit/repos/yacman/yacman/yacman.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">yacman.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/pepkit/repos/yacman/yacman/yacman.py#674" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">674</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> No output schema supplied. Running in schema-optional mode <span style="font-weight: bold">(</span><span style="color: #808000; text-decoration-color: #808000">validate_results</span>=<span style="color: #ff0000; text-decoration-color: #ff0000; font-style: italic">False</span><span style="font-weight: bold">)</span>.      <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/pipestat/pipestat.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">pipestat.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/pipestat/pipestat.py#1008" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1008</span></a>
         Results will not be validated against a schema.                                           <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">                </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Initializing results file                                                               <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/pipestat/backends/file_backend/filebackend.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">filebackend.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/pipestat/backends/file_backend/filebackend.py#785" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">785</span></a>
         <span style="color: #008000; text-decoration-color: #008000">'/tmp/refgenie_demo_2hzk45es/builds/t7/bowtie2_index/2.3.5/stats.yaml'</span>                  <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">                  </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> File does not exist, but create_file is true. Creating<span style="color: #808000; text-decoration-color: #808000">...</span>                              <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/ubiquerg/file_locking.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">file_locking.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/ubiquerg/file_locking.py#313" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">313</span></a>
</pre>



    ### Pipeline run code and environment:
    


    *          Command: `/home/nsheff/.local/lib/python3.12/site-packages/ipykernel_launcher.py -f /tmp/tmp30pcozza.json --HistoryManager.hist_file=:memory:`


    *     Compute host: `zither`


    *      Working dir: `/home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/docs`


    *        Outfolder: `/tmp/refgenie_demo_2hzk45es/builds/t7/bowtie2_index/2.3.5/`


    *         Log file: `/tmp/refgenie_demo_2hzk45es/builds/t7/bowtie2_index/2.3.5/refgenie_t7_bowtie2_index_2.3.5_log.md`


    *       Start time:  (08-20 21:55:45) elapsed: 0.0 _TIME_


    
    ### Version log:
    


    *   Python version: `3.12.3`


    *      Pypiper dir: `/home/nsheff/.local/lib/python3.12/site-packages/pypiper`


    *  Pypiper version: `0.15.1`


    * Pipestat version: `0.13.1`


    *     Pipeline dir: `/home/nsheff/.local/lib/python3.12/site-packages`


    * Pipeline version: 


    
    ### Arguments passed to pipeline:
    


    
    ### Initialized Pipestat Object:
    


    * PipestatManager (refgenie_t7_bowtie2_index_2.3.5)


    * Backend: File


    *  - results: /tmp/refgenie_demo_2hzk45es/builds/t7/bowtie2_index/2.3.5/stats.yaml


    *  - status: /tmp/refgenie_demo_2hzk45es/builds/t7/bowtie2_index/2.3.5


    * Multiple Pipelines Allowed: False


    * Pipeline name: refgenie_t7_bowtie2_index_2.3.5


    * Pipeline type: sample


    * Status Schema key: None


    * Results formatter: default_formatter


    * Results schema source: None


    * Status schema source: None


    * Records count: 0


    * Sample name: DEFAULT_SAMPLE_NAME
    


    
    ----------------------------------------
    


    Target to produce: `/tmp/refgenie_demo_2hzk45es/builds/t7/bowtie2_index/2.3.5/t7_bowtie2_index__2.3.5.flag`  


    
    > `bowtie2-build --threads 8 /tmp/refgenie_demo_2hzk45es/data/kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB/fasta/949c705a1c1668f0ea3d7da883d1a65e902b945b67de576407349ff72a8b1bf1/kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB.fa /tmp/refgenie_demo_2hzk45es/data/kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB/bowtie2_index/2.3.5/kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB` (1045591)


    <pre>


    [90m[[0m2026-08-21T01:55:45Z [33mWARN [0m bulker::shimlink[90m][0m 'bowtie2-build': `docker_command` is deprecated; use the `entrypoint` manifest field instead


    Settings:


    Building a SMALL index
      Output files: "/tmp/refgenie_demo_2hzk45es/data/kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB/bowtie2_index/2.3.5/kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB.*.bt2"


      Line rate: 6 (line is 64 bytes)


      Lines per side: 1 (side is 64 bytes)


      Offset rate: 4 (one in 16)


      FTable chars: 10


      Strings: unpacked


      Max bucket size: default


      Max bucket size, sqrt multiplier: default


      Max bucket size, len divisor: 32


      Difference-cover sample period: 1024


      Endianness: little


      Actual local endianness: little


      Sanity checking: disabled


      Assertions: disabled


      Random seed: 0


      Sizeofs: void*:8, int:4, long:8, size_t:8


    Input files DNA, FASTA:


      /tmp/refgenie_demo_2hzk45es/data/kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB/fasta/949c705a1c1668f0ea3d7da883d1a65e902b945b67de576407349ff72a8b1bf1/kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB.fa


    Reading reference sizes


      Time reading reference sizes: 00:00:00


    Calculating joined length


    Writing header


    Reserving space for joined string


    Joining reference sequences


      Time to join reference sequences: 00:00:00


    bmax according to bmaxDivN setting: 1248


    Using parameters --bmax 936 --dcv 1024


      Doing ahead-of-time memory usage test


      Passed!  Constructing with these parameters: --bmax 936 --dcv 1024


    Constructing suffix-array element generator


    Building DifferenceCoverSample


      Building sPrime


      Building sPrimeOrder


      V-Sorting samples


      V-Sorting samples time: 00:00:00


      Allocating rank array


      Ranking v-sort output


      Ranking v-sort output time: 00:00:00


      Invoking Larsson-Sadakane on ranks


      Invoking Larsson-Sadakane on ranks time: 00:00:00


      Sanity-checking and returning


    Building samples


    Reserving space for 86 sample suffixes


    Generating random suffixes


    QSorting 86 sample offsets, eliminating duplicates


    QSorting sample offsets, eliminating duplicates time: 00:00:00


    Multikey QSorting 86 samples


      (Using difference cover)


      Multikey QSorting samples time: 00:00:00


    Calculating bucket sizes


    Splitting and merging


      Splitting and merging time: 00:00:00


    Split 13, merged 35; iterating...


    Splitting and merging


      Splitting and merging time: 00:00:00


    Split 7, merged 7; iterating...


    Splitting and merging


      Splitting and merging time: 00:00:00


    Split 2, merged 4; iterating...


    Splitting and merging


      Splitting and merging time: 00:00:00


    Avg bucket size: 653.721 (target: 935)


    Converting suffix-array elements to index image


    Allocating ftab, absorbFtab


    Entering Ebwt loop


    Getting block 1 of 61


      Reserving size (936) for bucket 1


    Getting block 2 of 61


    Getting block 3 of 61


    Getting block 4 of 61


    Getting block 5 of 61


    Getting block 6 of 61


    Getting block 7 of 61


    Getting block 8 of 61


      Calculating Z arrays for bucket 1


      Reserving size (936) for bucket 2


      Reserving size (936) for bucket 3


      Reserving size (936) for bucket 4


      Reserving size (936) for bucket 5


      Reserving size (936) for bucket 6


      Reserving size (936) for bucket 7


      Reserving size (936) for bucket 8


      Entering block accumulator loop for bucket 1:


      Calculating Z arrays for bucket 2


      Calculating Z arrays for bucket 3


      Calculating Z arrays for bucket 4


      Calculating Z arrays for bucket 5


      Calculating Z arrays for bucket 6


      Calculating Z arrays for bucket 7


      Calculating Z arrays for bucket 8


      Entering block accumulator loop for bucket 2:


      Entering block accumulator loop for bucket 3:


      Entering block accumulator loop for bucket 5:


      Entering block accumulator loop for bucket 4:


      Entering block accumulator loop for bucket 8:


      Entering block accumulator loop for bucket 7:


      Entering block accumulator loop for bucket 6:


      bucket 1: 10%


      bucket 5: 10%


      bucket 1: 20%


      bucket 2: 10%


      bucket 4: 10%


      bucket 3: 10%


      bucket 8: 10%


      bucket 7: 10%


      bucket 6: 10%


      bucket 5: 20%


      bucket 1: 30%


      bucket 2: 20%


      bucket 4: 20%


      bucket 3: 20%


      bucket 7: 20%


      bucket 8: 20%


      bucket 5: 30%


      bucket 1: 40%


      bucket 6: 20%


      bucket 2: 30%


      bucket 4: 30%


      bucket 7: 30%


      bucket 3: 30%


      bucket 8: 30%


      bucket 5: 40%


      bucket 1: 50%


      bucket 6: 30%


      bucket 2: 40%


      bucket 4: 40%


      bucket 7: 40%


      bucket 3: 40%


      bucket 8: 40%


      bucket 1: 60%


      bucket 5: 50%


      bucket 6: 40%


      bucket 2: 50%


      bucket 4: 50%


      bucket 7: 50%


      bucket 1: 70%


      bucket 3: 50%


      bucket 5: 60%


      bucket 8: 50%


      bucket 6: 50%


      bucket 1: 80%


      bucket 2: 60%


      bucket 7: 60%


      bucket 4: 60%


      bucket 3: 60%


      bucket 5: 70%


      bucket 8: 60%


      bucket 1: 90%


      bucket 6: 60%


      bucket 2: 70%


      bucket 7: 70%


      bucket 4: 70%


      bucket 3: 70%


      bucket 5: 80%


      bucket 1: 100%


      bucket 8: 70%


      bucket 6: 70%


      bucket 7: 80%


      bucket 2: 80%


      bucket 4: 80%


      Sorting block of length 878 for bucket 1


      (Using difference cover)


      bucket 3: 80%


      bucket 5: 90%


      bucket 8: 80%


      bucket 6: 80%


      bucket 7: 90%


      bucket 2: 90%


      bucket 4: 90%


      bucket 3: 90%


      bucket 5: 100%


      bucket 8: 90%


      bucket 7: 100%


      bucket 6: 90%


      bucket 2: 100%


      Sorting block of length 395 for bucket 5


      (Using difference cover)


      bucket 4: 100%


      Sorting block time: 00:00:00


      Sorting block time: 00:00:00


      Sorting block of length 709 for bucket 7


      (Using difference cover)


      bucket 3: 100%


      Sorting block of length 639 for bucket 2


      (Using difference cover)


      bucket 8: 100%


      Sorting block of length 552 for bucket 4


      (Using difference cover)


    Returning block of 879 for bucket 1


      bucket 6: 100%


    Returning block of 396 for bucket 5


      Sorting block of length 516 for bucket 3


      (Using difference cover)


      Sorting block of length 598 for bucket 8


      (Using difference cover)


      Sorting block of length 596 for bucket 6


      (Using difference cover)


    Getting block 9 of 61


      Sorting block time: 00:00:00


    Getting block 10 of 61


      Reserving size (936) for bucket 9


    Returning block of 710 for bucket 7


      Reserving size (936) for bucket 10


      Calculating Z arrays for bucket 9


      Calculating Z arrays for bucket 10


      Entering block accumulator loop for bucket 9:


      Entering block accumulator loop for bucket 10:


    Getting block 11 of 61


      bucket 10: 10%


      bucket 9: 10%


      Reserving size (936) for bucket 11


      Calculating Z arrays for bucket 11


      Sorting block time: 00:00:00


    Returning block of 599 for bucket 8


      Sorting block time: 00:00:00


      Entering block accumulator loop for bucket 11:


    Returning block of 597 for bucket 6


      bucket 10: 20%


      bucket 9: 20%


    Getting block 12 of 61


      Reserving size (936) for bucket 12


      Calculating Z arrays for bucket 12


      Entering block accumulator loop for bucket 12:


      bucket 11: 10%


      Sorting block time: 00:00:00


      bucket 9: 30%


      bucket 10: 30%


      Sorting block time: 00:00:00


    Returning block of 640 for bucket 2


      Sorting block time: 00:00:00


    Getting block 13 of 61


      bucket 11: 20%


      bucket 12: 10%


      bucket 9: 40%


    Returning block of 553 for bucket 4


    Returning block of 517 for bucket 3


      Reserving size (936) for bucket 13


      bucket 10: 40%


      bucket 11: 30%


      Calculating Z arrays for bucket 13


      bucket 12: 20%


      bucket 9: 50%


    Getting block 14 of 61


      Entering block accumulator loop for bucket 13:


      bucket 10: 50%


      bucket 11: 40%


      bucket 9: 60%


      bucket 12: 30%


      Reserving size (936) for bucket 14


    Getting block 15 of 61


    Getting block 16 of 61


      bucket 10: 60%


      bucket 9: 70%


      bucket 11: 50%


      bucket 12: 40%


      Calculating Z arrays for bucket 14


      Reserving size (936) for bucket 15


      Reserving size (936) for bucket 16


      bucket 13: 10%


      bucket 9: 80%


      bucket 10: 70%


      bucket 11: 60%


      bucket 12: 50%


      Entering block accumulator loop for bucket 14:


      Calculating Z arrays for bucket 15


      Calculating Z arrays for bucket 16


      bucket 9: 90%


      bucket 13: 20%


      bucket 10: 80%


      bucket 11: 70%


      bucket 12: 60%


      Entering block accumulator loop for bucket 15:


      Entering block accumulator loop for bucket 16:


      bucket 9: 100%


      bucket 12: 70%


      bucket 13: 30%


      bucket 11: 80%


      bucket 10: 90%


      Sorting block of length 581 for bucket 9


      (Using difference cover)


      bucket 14: 10%


      Sorting block time: 00:00:00


      bucket 12: 80%


      bucket 11: 90%


      bucket 13: 40%


      bucket 10: 100%


      bucket 16: 10%


      bucket 15: 10%


    Returning block of 582 for bucket 9


      bucket 14: 20%


      Sorting block of length 703 for bucket 10


      (Using difference cover)


      bucket 11: 100%


      bucket 12: 90%


      bucket 13: 50%


      bucket 16: 20%


      bucket 15: 20%


      Sorting block time: 00:00:00


      bucket 14: 30%


      Sorting block of length 457 for bucket 11


      (Using difference cover)


      bucket 12: 100%


      bucket 13: 60%


      Sorting block time: 00:00:00


    Getting block 17 of 61


    Returning block of 704 for bucket 10


      bucket 16: 30%


      bucket 15: 30%


      bucket 14: 40%


      Sorting block of length 686 for bucket 12


      (Using difference cover)


    Returning block of 458 for bucket 11


      Reserving size (936) for bucket 17


      bucket 13: 70%


      bucket 16: 40%


      bucket 15: 40%


      Sorting block time: 00:00:00


      bucket 14: 50%


      Calculating Z arrays for bucket 17


    Getting block 18 of 61


      bucket 13: 80%


      bucket 16: 50%


      bucket 15: 50%


    Returning block of 687 for bucket 12


      Reserving size (936) for bucket 18


    Getting block 19 of 61


      Entering block accumulator loop for bucket 17:


      bucket 14: 60%


      Calculating Z arrays for bucket 18


      Reserving size (936) for bucket 19


      bucket 16: 60%


      bucket 13: 90%


      bucket 14: 70%


      bucket 15: 60%


      Entering block accumulator loop for bucket 18:


      Calculating Z arrays for bucket 19


    Getting block 20 of 61


      bucket 17: 10%


      Entering block accumulator loop for bucket 19:


      bucket 16: 70%


      bucket 13: 100%


      bucket 14: 80%


      bucket 15: 70%


      Reserving size (936) for bucket 20


      Sorting block of length 595 for bucket 13


      (Using difference cover)


      bucket 17: 20%


      bucket 16: 80%


      bucket 18: 10%


      bucket 14: 90%


      bucket 15: 80%


      Sorting block time: 00:00:00


      Calculating Z arrays for bucket 20


      bucket 19: 10%


      bucket 17: 30%


      bucket 16: 90%


    Returning block of 596 for bucket 13


      bucket 14: 100%


      bucket 18: 20%


      bucket 15: 90%


      Entering block accumulator loop for bucket 20:


      Sorting block of length 402 for bucket 14


      (Using difference cover)


      bucket 19: 20%


      bucket 17: 40%


      bucket 16: 100%


      Sorting block time: 00:00:00


      Sorting block of length 768 for bucket 16


      (Using difference cover)


      bucket 18: 30%


      bucket 15: 100%


      bucket 17: 50%


      bucket 19: 30%


    Returning block of 403 for bucket 14


    Getting block 21 of 61


      Sorting block time: 00:00:00


      bucket 20: 10%


      Sorting block of length 819 for bucket 15


      (Using difference cover)


      bucket 18: 40%


      bucket 17: 60%


      bucket 19: 40%


      Reserving size (936) for bucket 21


    Returning block of 769 for bucket 16


      bucket 20: 20%


      Sorting block time: 00:00:00


      bucket 18: 50%


      Calculating Z arrays for bucket 21


      bucket 17: 70%


      bucket 19: 50%


    Getting block 22 of 61


    Returning block of 820 for bucket 15


      Entering block accumulator loop for bucket 21:


      bucket 20: 30%


      bucket 17: 80%


      bucket 18: 60%


      bucket 19: 60%


      Reserving size (936) for bucket 22


    Getting block 23 of 61


      bucket 20: 40%


      bucket 17: 90%


      Calculating Z arrays for bucket 22


    Getting block 24 of 61


      Reserving size (936) for bucket 23


      bucket 19: 70%


      Entering block accumulator loop for bucket 22:


      bucket 21: 10%


      bucket 20: 50%


      bucket 17: 100%


      Reserving size (936) for bucket 24


      Calculating Z arrays for bucket 23


      bucket 19: 80%


      Sorting block of length 736 for bucket 17


      (Using difference cover)


      bucket 21: 20%


      Calculating Z arrays for bucket 24


      bucket 20: 60%


      bucket 18: 70%


      Sorting block time: 00:00:00


      Entering block accumulator loop for bucket 23:


      bucket 22: 10%


      bucket 19: 90%


      Entering block accumulator loop for bucket 24:


      bucket 21: 30%


      bucket 20: 70%


    Returning block of 737 for bucket 17


      bucket 22: 20%


      bucket 19: 100%


      bucket 21: 40%


      bucket 20: 80%


      bucket 18: 80%


      Sorting block of length 902 for bucket 19


      (Using difference cover)


      bucket 24: 10%


      bucket 22: 30%


      bucket 21: 50%


      bucket 20: 90%


      bucket 23: 10%


      Sorting block time: 00:00:00


      bucket 18: 90%


    Getting block 25 of 61


      bucket 24: 20%


      bucket 22: 40%


      bucket 21: 60%


      bucket 20: 100%


    Returning block of 903 for bucket 19


      bucket 23: 20%


      Reserving size (936) for bucket 25


      bucket 18: 100%


      bucket 24: 30%


      bucket 22: 50%


      Sorting block of length 578 for bucket 20


      (Using difference cover)


      bucket 21: 70%


      Calculating Z arrays for bucket 25


      bucket 23: 30%


      Sorting block time: 00:00:00


      Sorting block of length 693 for bucket 18


      (Using difference cover)


      bucket 24: 40%


      bucket 22: 60%


    Getting block 26 of 61


      Entering block accumulator loop for bucket 25:


      bucket 21: 80%


      Sorting block time: 00:00:00


    Returning block of 579 for bucket 20


      bucket 23: 40%


      bucket 24: 50%


      Reserving size (936) for bucket 26


      bucket 22: 70%


    Returning block of 694 for bucket 18


      bucket 21: 90%


      Calculating Z arrays for bucket 26


      bucket 24: 60%


      bucket 25: 10%


      bucket 22: 80%


    Getting block 27 of 61


      bucket 23: 50%


      Entering block accumulator loop for bucket 26:


      bucket 21: 100%


      bucket 24: 70%


      Reserving size (936) for bucket 27


      bucket 25: 20%


      bucket 22: 90%


      Sorting block of length 599 for bucket 21


      (Using difference cover)


      bucket 23: 60%


      Calculating Z arrays for bucket 27


      bucket 24: 80%


    Getting block 28 of 61


      Sorting block time: 00:00:00


      bucket 25: 30%


      bucket 22: 100%


      bucket 26: 10%


      Entering block accumulator loop for bucket 27:


      bucket 23: 70%


      Reserving size (936) for bucket 28


    Returning block of 600 for bucket 21


      bucket 24: 90%


      Sorting block of length 895 for bucket 22


      (Using difference cover)


      bucket 25: 40%


      bucket 26: 20%


      Calculating Z arrays for bucket 28


      bucket 23: 80%


      bucket 25: 50%


      bucket 24: 100%


      bucket 26: 30%


      Entering block accumulator loop for bucket 28:


      bucket 27: 10%


    Getting block 29 of 61


      bucket 23: 90%


      Sorting block of length 583 for bucket 24


      (Using difference cover)


      Sorting block time: 00:00:00


      Reserving size (936) for bucket 29


    Returning block of 896 for bucket 22


      bucket 25: 60%


      bucket 26: 40%


      bucket 27: 20%


      bucket 23: 100%


      bucket 28: 10%


      Sorting block time: 00:00:00


      Calculating Z arrays for bucket 29


      Sorting block of length 917 for bucket 23


      (Using difference cover)


    Returning block of 584 for bucket 24


      bucket 25: 70%


      bucket 27: 30%


      bucket 26: 50%


      bucket 28: 20%


      Entering block accumulator loop for bucket 29:


    Getting block 30 of 61


      bucket 25: 80%


      bucket 26: 60%


      bucket 27: 40%


      bucket 28: 30%


      Sorting block time: 00:00:00


      Reserving size (936) for bucket 30


    Getting block 31 of 61


      bucket 25: 90%


    Returning block of 918 for bucket 23


      bucket 26: 70%


      bucket 27: 50%


      bucket 28: 40%


      Calculating Z arrays for bucket 30


      Reserving size (936) for bucket 31


      bucket 29: 10%


      bucket 25: 100%


      Calculating Z arrays for bucket 31


      Entering block accumulator loop for bucket 30:


      bucket 26: 80%


      bucket 27: 60%


      bucket 28: 50%


      Sorting block of length 490 for bucket 25


      (Using difference cover)


    Getting block 32 of 61


      Entering block accumulator loop for bucket 31:


      bucket 29: 20%


      Reserving size (936) for bucket 32


      bucket 26: 90%


      bucket 27: 70%


      Sorting block time: 00:00:00


      bucket 28: 60%


      bucket 30: 10%


      Calculating Z arrays for bucket 32


    Returning block of 491 for bucket 25


      bucket 29: 30%


      bucket 26: 100%


      bucket 31: 10%


      bucket 27: 80%


      Entering block accumulator loop for bucket 32:


      bucket 28: 70%


      bucket 30: 20%


      Sorting block of length 625 for bucket 26


      (Using difference cover)


      bucket 29: 40%


    Getting block 33 of 61


      bucket 27: 90%


      bucket 31: 20%


      Sorting block time: 00:00:00


      bucket 28: 80%


      Reserving size (936) for bucket 33


      bucket 30: 30%


    Returning block of 626 for bucket 26


      bucket 32: 10%


      bucket 29: 50%


      Calculating Z arrays for bucket 33


      bucket 27: 100%


      bucket 31: 30%


      bucket 28: 90%


      bucket 30: 40%


      Sorting block of length 532 for bucket 27


      (Using difference cover)


    Getting block 34 of 61


      bucket 32: 20%


      Entering block accumulator loop for bucket 33:


      bucket 29: 60%


      bucket 31: 40%


      Reserving size (936) for bucket 34


      Sorting block time: 00:00:00


      bucket 28: 100%


      bucket 30: 50%


      Calculating Z arrays for bucket 34


    Returning block of 533 for bucket 27


      bucket 32: 30%


      bucket 29: 70%


      Sorting block of length 775 for bucket 28


      (Using difference cover)


      bucket 31: 50%


      bucket 33: 10%


      Entering block accumulator loop for bucket 34:


      bucket 30: 60%


      bucket 32: 40%


      Sorting block time: 00:00:00


      bucket 29: 80%


      bucket 31: 60%


    Getting block 35 of 61


      bucket 33: 20%


    Returning block of 776 for bucket 28


      bucket 30: 70%


      Reserving size (936) for bucket 35


      bucket 32: 50%


      bucket 31: 70%


      bucket 33: 30%


      bucket 29: 90%


      bucket 34: 10%


      Calculating Z arrays for bucket 35


      bucket 30: 80%


      bucket 32: 60%


      Entering block accumulator loop for bucket 35:


      bucket 31: 80%


      bucket 33: 40%


      bucket 29: 100%


    Getting block 36 of 61


      bucket 34: 20%


      bucket 30: 90%


      bucket 32: 70%


      Sorting block of length 397 for bucket 29


      (Using difference cover)


      Sorting block time: 00:00:00


      Reserving size (936) for bucket 36


      bucket 31: 90%


      bucket 33: 50%


      bucket 34: 30%


      bucket 32: 80%


      bucket 35: 10%


      bucket 30: 100%


    Returning block of 398 for bucket 29


      Calculating Z arrays for bucket 36


      bucket 31: 100%


      bucket 33: 60%


      bucket 34: 40%


      Sorting block of length 837 for bucket 30


      (Using difference cover)


      bucket 32: 90%


      bucket 35: 20%


      Sorting block of length 852 for bucket 31


      (Using difference cover)


      Entering block accumulator loop for bucket 36:


      bucket 33: 70%


      bucket 34: 50%


    Getting block 37 of 61


      Sorting block time: 00:00:00


      bucket 32: 100%


      Sorting block time: 00:00:00


      bucket 35: 30%


      Reserving size (936) for bucket 37


    Returning block of 838 for bucket 30


      bucket 33: 80%


      Sorting block of length 542 for bucket 32


      (Using difference cover)


      bucket 34: 60%


    Returning block of 853 for bucket 31


      Calculating Z arrays for bucket 37


      Sorting block time: 00:00:00


      bucket 36: 10%


      bucket 35: 40%


      bucket 33: 90%


      bucket 34: 70%


      Entering block accumulator loop for bucket 37:


    Returning block of 543 for bucket 32


    Getting block 38 of 61


      bucket 36: 20%


      bucket 35: 50%


      bucket 33: 100%


    Getting block 39 of 61


      bucket 34: 80%


      Reserving size (936) for bucket 38


      Sorting block of length 571 for bucket 33


      (Using difference cover)


      bucket 37: 10%


      bucket 36: 30%


      bucket 35: 60%


      Sorting block time: 00:00:00


      Reserving size (936) for bucket 39


      Calculating Z arrays for bucket 38


      bucket 34: 90%


    Getting block 40 of 61


      bucket 37: 20%


      bucket 36: 40%


      bucket 35: 70%


    Returning block of 572 for bucket 33


      Calculating Z arrays for bucket 39


      Entering block accumulator loop for bucket 38:


      Reserving size (936) for bucket 40


      bucket 34: 100%


      bucket 37: 30%


      bucket 36: 50%


      bucket 35: 80%


      Entering block accumulator loop for bucket 39:


      Calculating Z arrays for bucket 40


      Sorting block of length 581 for bucket 34


      (Using difference cover)


    Getting block 41 of 61


      bucket 37: 40%


      bucket 38: 10%


      bucket 36: 60%


      bucket 35: 90%


      Sorting block time: 00:00:00


      Entering block accumulator loop for bucket 40:


      Reserving size (936) for bucket 41


      bucket 37: 50%


    Returning block of 582 for bucket 34


      bucket 36: 70%


      bucket 38: 20%


      Calculating Z arrays for bucket 41


      bucket 35: 100%


      bucket 39: 10%


      Sorting block of length 593 for bucket 35


      (Using difference cover)


      Entering block accumulator loop for bucket 41:


      bucket 37: 60%


      bucket 36: 80%


      bucket 38: 30%


      bucket 40: 10%


    Getting block 42 of 61


      Sorting block time: 00:00:00


      bucket 39: 20%


      bucket 37: 70%


      Reserving size (936) for bucket 42


    Returning block of 594 for bucket 35


      bucket 36: 90%


      bucket 40: 20%


      bucket 38: 40%


      bucket 41: 10%


      Calculating Z arrays for bucket 42


      bucket 39: 30%


      bucket 37: 80%


      Entering block accumulator loop for bucket 42:


      bucket 36: 100%


      bucket 40: 30%


      bucket 38: 50%


    Getting block 43 of 61


      bucket 41: 20%


      bucket 37: 90%


      bucket 39: 40%


      Sorting block of length 379 for bucket 36


      (Using difference cover)


      Reserving size (936) for bucket 43


      bucket 40: 40%


      bucket 42: 10%


      bucket 38: 60%


      bucket 41: 30%


      Sorting block time: 00:00:00


      bucket 37: 100%


      bucket 39: 50%


      Calculating Z arrays for bucket 43


    Returning block of 380 for bucket 36


      bucket 42: 20%


      bucket 38: 70%


      bucket 41: 40%


      Sorting block of length 600 for bucket 37


      (Using difference cover)


      Entering block accumulator loop for bucket 43:


      bucket 39: 60%


      bucket 42: 30%


      bucket 38: 80%


      bucket 39: 70%


      Sorting block time: 00:00:00


      bucket 41: 50%


    Getting block 44 of 61


    Returning block of 601 for bucket 37


      bucket 43: 10%


      bucket 42: 40%


      bucket 38: 90%


      bucket 39: 80%


      bucket 41: 60%


      Reserving size (936) for bucket 44


      bucket 43: 20%


      bucket 42: 50%


      bucket 38: 100%


      bucket 41: 70%


      bucket 39: 90%


      Calculating Z arrays for bucket 44


    Getting block 45 of 61


      Sorting block of length 697 for bucket 38


      (Using difference cover)


      bucket 43: 30%


      bucket 42: 60%


      bucket 41: 80%


      bucket 39: 100%


      Reserving size (936) for bucket 45


      Sorting block time: 00:00:00


      Entering block accumulator loop for bucket 44:


      bucket 43: 40%


      Sorting block of length 889 for bucket 39


      (Using difference cover)


      bucket 42: 70%


      Calculating Z arrays for bucket 45


      bucket 41: 90%


    Returning block of 698 for bucket 38


      Entering block accumulator loop for bucket 45:


      bucket 43: 50%


      bucket 42: 80%


      Sorting block time: 00:00:00


      bucket 41: 100%


      bucket 44: 10%


    Returning block of 890 for bucket 39


      Sorting block of length 892 for bucket 41


      (Using difference cover)


    Getting block 46 of 61


      bucket 43: 60%


      bucket 42: 90%


      bucket 45: 10%


      Sorting block time: 00:00:00


      bucket 44: 20%


      Reserving size (936) for bucket 46


      bucket 43: 70%


      bucket 42: 100%


      bucket 45: 20%


    Returning block of 893 for bucket 41


      Calculating Z arrays for bucket 46


    Getting block 47 of 61


      Sorting block of length 652 for bucket 42


      (Using difference cover)


      bucket 44: 30%


      bucket 43: 80%


      Reserving size (936) for bucket 47


      bucket 45: 30%


      Entering block accumulator loop for bucket 46:


      Sorting block time: 00:00:00


      Calculating Z arrays for bucket 47


    Returning block of 653 for bucket 42


      bucket 44: 40%


    Getting block 48 of 61


      bucket 43: 90%


      bucket 45: 40%


      Entering block accumulator loop for bucket 47:


      bucket 46: 10%


      bucket 44: 50%


      Reserving size (936) for bucket 48


    Getting block 49 of 61


      bucket 43: 100%


      bucket 45: 50%


      bucket 46: 20%


      bucket 44: 60%


      Calculating Z arrays for bucket 48


      Reserving size (936) for bucket 49


      Sorting block of length 811 for bucket 43


      (Using difference cover)


      bucket 47: 10%


      bucket 45: 60%


      bucket 46: 30%


      bucket 44: 70%


      Calculating Z arrays for bucket 49


      Entering block accumulator loop for bucket 48:


      Sorting block time: 00:00:00


      Entering block accumulator loop for bucket 49:


      bucket 45: 70%


      bucket 47: 20%


      bucket 46: 40%


      bucket 44: 80%


    Returning block of 812 for bucket 43


      bucket 45: 80%


      bucket 47: 30%


      bucket 46: 50%


      bucket 44: 90%


      bucket 48: 10%


      bucket 49: 10%


      bucket 47: 40%


      bucket 45: 90%


      bucket 46: 60%


      bucket 44: 100%


    Getting block 50 of 61


      bucket 40: 50%


      bucket 48: 20%


      Sorting block of length 463 for bucket 44


      (Using difference cover)


      bucket 49: 20%


      Reserving size (936) for bucket 50


      bucket 47: 50%


      bucket 45: 100%


      bucket 46: 70%


      bucket 40: 60%


      Sorting block time: 00:00:00


      bucket 48: 30%


      Calculating Z arrays for bucket 50


      Sorting block of length 749 for bucket 45


      (Using difference cover)


      bucket 49: 30%


      bucket 47: 60%


      bucket 46: 80%


    Returning block of 464 for bucket 44


      bucket 40: 70%


      Entering block accumulator loop for bucket 50:


      bucket 48: 40%


      bucket 49: 40%


      bucket 47: 70%


      Sorting block time: 00:00:00


      bucket 46: 90%


      bucket 40: 80%


    Returning block of 750 for bucket 45


      bucket 48: 50%


    Getting block 51 of 61


      bucket 49: 50%


      bucket 50: 10%


      bucket 47: 80%


      bucket 40: 90%


      bucket 48: 60%


      Reserving size (936) for bucket 51


      bucket 49: 60%


      bucket 40: 100%


      bucket 48: 70%


      bucket 50: 20%


      Calculating Z arrays for bucket 51


      bucket 47: 90%


      Sorting block of length 619 for bucket 40


      (Using difference cover)


      bucket 49: 70%


      bucket 48: 80%


      bucket 50: 30%


      Entering block accumulator loop for bucket 51:


      bucket 49: 80%


      bucket 48: 90%


      Sorting block time: 00:00:00


      bucket 47: 100%


      bucket 50: 40%


      bucket 49: 90%


    Returning block of 620 for bucket 40


      bucket 48: 100%


      Sorting block of length 928 for bucket 47


      (Using difference cover)


      bucket 51: 10%


      bucket 50: 50%


      bucket 46: 100%


      Sorting block of length 419 for bucket 48


      (Using difference cover)


      bucket 49: 100%


      Sorting block of length 368 for bucket 46


      (Using difference cover)


      bucket 50: 60%


      Sorting block time: 00:00:00


      Sorting block time: 00:00:00


      bucket 51: 20%


      Sorting block of length 893 for bucket 49


      (Using difference cover)


    Getting block 52 of 61


      Sorting block time: 00:00:00


    Getting block 53 of 61


      Sorting block time: 00:00:00


    Returning block of 929 for bucket 47


    Returning block of 420 for bucket 48


      bucket 50: 70%


      bucket 51: 30%


      Reserving size (936) for bucket 52


    Returning block of 369 for bucket 46


    Returning block of 894 for bucket 49


      Reserving size (936) for bucket 53


      bucket 50: 80%


      Calculating Z arrays for bucket 52


      Calculating Z arrays for bucket 53


      bucket 51: 40%


    Getting block 54 of 61


    Getting block 55 of 61


      Entering block accumulator loop for bucket 52:


    Getting block 56 of 61


      bucket 50: 90%


      Entering block accumulator loop for bucket 53:


      bucket 51: 50%


      Reserving size (936) for bucket 54


      Reserving size (936) for bucket 55


    Getting block 57 of 61


      Reserving size (936) for bucket 56


      bucket 50: 100%


      bucket 51: 60%


      Calculating Z arrays for bucket 54


      Calculating Z arrays for bucket 55


      Reserving size (936) for bucket 57


      Calculating Z arrays for bucket 56


      Sorting block of length 218 for bucket 50


      (Using difference cover)


      bucket 52: 10%


      bucket 53: 10%


      bucket 51: 70%


      Calculating Z arrays for bucket 57


      Entering block accumulator loop for bucket 54:


      Sorting block time: 00:00:00


      Entering block accumulator loop for bucket 55:


      Entering block accumulator loop for bucket 56:


      Entering block accumulator loop for bucket 57:


      bucket 52: 20%


      bucket 51: 80%


      bucket 53: 20%


    Returning block of 219 for bucket 50


      bucket 52: 30%


      bucket 51: 90%


      bucket 53: 30%


      bucket 55: 10%


      bucket 54: 10%


      bucket 56: 10%


      bucket 57: 10%


      bucket 52: 40%


    Getting block 58 of 61


      bucket 51: 100%


      bucket 53: 40%


      bucket 55: 20%


      Reserving size (936) for bucket 58


      bucket 54: 20%


      Sorting block of length 889 for bucket 51


      (Using difference cover)


      bucket 56: 20%


      Calculating Z arrays for bucket 58


      bucket 57: 20%


      bucket 52: 50%


      bucket 55: 30%


      bucket 53: 50%


      Sorting block time: 00:00:00


      bucket 54: 30%


      Entering block accumulator loop for bucket 58:


      bucket 56: 30%


      bucket 57: 30%


      bucket 52: 60%


      bucket 55: 40%


    Returning block of 890 for bucket 51


      bucket 53: 60%


      bucket 54: 40%


      bucket 56: 40%


      bucket 57: 40%


      bucket 52: 70%


      bucket 55: 50%


      bucket 58: 10%


      bucket 53: 70%


      bucket 54: 50%


      bucket 56: 50%


      bucket 57: 50%


    Getting block 59 of 61


      bucket 52: 80%


      bucket 55: 60%


      bucket 53: 80%


      bucket 58: 20%


      Reserving size (936) for bucket 59


      bucket 57: 60%


      bucket 54: 60%


      bucket 56: 60%


      bucket 52: 90%


      bucket 53: 90%


      bucket 55: 70%


      Calculating Z arrays for bucket 59


      bucket 58: 30%


      bucket 57: 70%


      bucket 53: 100%


      bucket 54: 70%


      bucket 55: 80%


      bucket 56: 70%


      bucket 52: 100%


      Entering block accumulator loop for bucket 59:


      Sorting block of length 934 for bucket 53


      (Using difference cover)


      Sorting block of length 845 for bucket 52


      Sorting block time: 00:00:00


      bucket 58: 40%


      bucket 57: 80%


      bucket 54: 80%


      bucket 55: 90%


      bucket 56: 80%


    Returning block of 935 for bucket 53


      bucket 59: 10%


      bucket 58: 50%


      bucket 57: 90%


      bucket 54: 90%


      bucket 55: 100%


      bucket 56: 90%


      Sorting block of length 897 for bucket 55


      (Using difference cover)


      bucket 59: 20%


      bucket 58: 60%


      bucket 57: 100%


    Getting block 60 of 61


      bucket 54: 100%


      Sorting block of length 493 for bucket 57


      (Using difference cover)


      bucket 56: 100%


      bucket 59: 30%


      Reserving size (936) for bucket 60


      Sorting block of length 62 for bucket 54


      (Using difference cover)


      bucket 58: 70%


      Sorting block of length 606 for bucket 56


      Sorting block time: 00:00:00


      (Using difference cover)


      bucket 59: 40%


      Sorting block time: 00:00:00


      Sorting block time: 00:00:00


      Calculating Z arrays for bucket 60


    Returning block of 494 for bucket 57


      bucket 58: 80%


    Returning block of 898 for bucket 55


    Returning block of 63 for bucket 54


      bucket 59: 50%


      Entering block accumulator loop for bucket 60:


      bucket 58: 90%


      bucket 59: 60%


      Sorting block time: 00:00:00


    Getting block 61 of 61


    Returning block of 607 for bucket 56


      bucket 59: 70%


      bucket 58: 100%


      Reserving size (936) for bucket 61


      bucket 60: 10%


      Sorting block of length 924 for bucket 58


      (Using difference cover)


      Calculating Z arrays for bucket 61


      bucket 59: 80%


      Entering block accumulator loop for bucket 61:


      bucket 60: 20%


      Sorting block time: 00:00:00


      bucket 59: 90%


    Returning block of 925 for bucket 58


      bucket 61: 10%


      bucket 60: 30%


      bucket 61: 20%


      bucket 59: 100%


      bucket 61: 30%


      bucket 60: 40%


      Sorting block of length 638 for bucket 59


      (Using difference cover)


      bucket 61: 40%


      Sorting block time: 00:00:00


      bucket 60: 50%


      bucket 61: 50%


    Returning block of 639 for bucket 59


      bucket 61: 60%


      (Using difference cover)


      bucket 60: 60%


      bucket 61: 70%


      bucket 61: 80%


      bucket 60: 70%


      Sorting block time: 00:00:00


    Returning block of 846 for bucket 52


      bucket 61: 90%


      bucket 60: 80%


      bucket 61: 100%


      Sorting block of length 693 for bucket 61


      (Using difference cover)


      bucket 60: 90%


      bucket 60: 100%


      Sorting block of length 727 for bucket 60


      (Using difference cover)


      Sorting block time: 00:00:00


    Returning block of 694 for bucket 61


      Sorting block time: 00:00:00


    Returning block of 728 for bucket 60


    Exited Ebwt loop


    fchr[A]: 0


    fchr[C]: 10842


    fchr[G]: 19880


    fchr[T]: 30171


    fchr[$]: 39937


    Exiting Ebwt::buildToDisk()


    Returning from initFromVector


    Wrote 4207850 bytes to primary EBWT file: /tmp/refgenie_demo_2hzk45es/data/kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB/bowtie2_index/2.3.5/kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB.1.bt2


    Wrote 9992 bytes to secondary EBWT file: /tmp/refgenie_demo_2hzk45es/data/kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB/bowtie2_index/2.3.5/kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB.2.bt2


    Re-opening _in1 and _in2 as input streams


    Returning from Ebwt constructor


    Headers:


        len: 39937


        bwtLen: 39938


        sz: 9985


        bwtSz: 9985


        lineRate: 6


        offRate: 4


        offMask: 0xfffffff0


        ftabChars: 10


        eftabLen: 20


        eftabSz: 80


        ftabLen: 1048577


        ftabSz: 4194308


        offsLen: 2497


        offsSz: 9988


        lineSz: 64


        sideSz: 64


        sideBwtSz: 48


        sideBwtLen: 192


        numSides: 209


        numLines: 209


        ebwtTotLen: 13376


        ebwtTotSz: 13376


        color: 0


        reverse: 0


    Total time for call to driver() for forward index: 00:00:00


    Reading reference sizes


      Time reading reference sizes: 00:00:00


    Calculating joined length


    Writing header


    Reserving space for joined string


    Joining reference sequences


      Time to join reference sequences: 00:00:00


      Time to reverse reference sequence: 00:00:00


    bmax according to bmaxDivN setting: 1248


    Using parameters --bmax 936 --dcv 1024


      Doing ahead-of-time memory usage test


      Passed!  Constructing with these parameters: --bmax 936 --dcv 1024


    Constructing suffix-array element generator


    Building DifferenceCoverSample


      Building sPrime


      Building sPrimeOrder


      V-Sorting samples


      V-Sorting samples time: 00:00:00


      Allocating rank array


      Ranking v-sort output


      Ranking v-sort output time: 00:00:00


      Invoking Larsson-Sadakane on ranks


      Invoking Larsson-Sadakane on ranks time: 00:00:00


      Sanity-checking and returning


    Building samples


    Reserving space for 86 sample suffixes


    Generating random suffixes


    QSorting 86 sample offsets, eliminating duplicates


    QSorting sample offsets, eliminating duplicates time: 00:00:00


    Multikey QSorting 86 samples


      (Using difference cover)


      Multikey QSorting samples time: 00:00:00


    Calculating bucket sizes


    Splitting and merging


      Splitting and merging time: 00:00:00


    Split 13, merged 40; iterating...


    Splitting and merging


      Splitting and merging time: 00:00:00


    Split 7, merged 6; iterating...


    Splitting and merging


      Splitting and merging time: 00:00:00


    Split 3, merged 4; iterating...


    Splitting and merging


      Splitting and merging time: 00:00:00


    Split 2, merged 3; iterating...


    Splitting and merging


      Splitting and merging time: 00:00:00


    Split 2, merged 1; iterating...


    Avg bucket size: 664.633 (target: 935)


    Converting suffix-array elements to index image


    Allocating ftab, absorbFtab


    Entering Ebwt loop


    Getting block 1 of 60


      Reserving size (936) for bucket 1


    Getting block 2 of 60


    Getting block 3 of 60


    Getting block 4 of 60


      Calculating Z arrays for bucket 1


    Getting block 5 of 60


      Reserving size (936) for bucket 2


    Getting block 6 of 60


      Reserving size (936) for bucket 3


      Reserving size (936) for bucket 4


    Getting block 7 of 60


      Entering block accumulator loop for bucket 1:


      Reserving size (936) for bucket 5


    Getting block 8 of 60


      Calculating Z arrays for bucket 2


      Reserving size (936) for bucket 6


      Calculating Z arrays for bucket 3


      Calculating Z arrays for bucket 4


      Reserving size (936) for bucket 7


      Calculating Z arrays for bucket 5


      Reserving size (936) for bucket 8


      Calculating Z arrays for bucket 6


      Entering block accumulator loop for bucket 3:


      Calculating Z arrays for bucket 7


      Entering block accumulator loop for bucket 2:


      Entering block accumulator loop for bucket 4:


      Calculating Z arrays for bucket 8


      Entering block accumulator loop for bucket 5:


      Entering block accumulator loop for bucket 6:


      bucket 1: 10%


      Entering block accumulator loop for bucket 7:


      Entering block accumulator loop for bucket 8:


      bucket 3: 10%


      bucket 4: 10%


      bucket 1: 20%


      bucket 5: 10%


      bucket 2: 10%


      bucket 6: 10%


      bucket 4: 20%


      bucket 3: 20%


      bucket 8: 10%


      bucket 1: 30%


      bucket 5: 20%


      bucket 7: 10%


      bucket 4: 30%


      bucket 2: 20%


      bucket 3: 30%


      bucket 1: 40%


      bucket 6: 20%


      bucket 8: 20%


      bucket 4: 40%


      bucket 5: 30%


      bucket 3: 40%


      bucket 1: 50%


      bucket 2: 30%


      bucket 6: 30%


      bucket 8: 30%


      bucket 4: 50%


      bucket 3: 50%


      bucket 7: 20%


      bucket 1: 60%


      bucket 5: 40%


      bucket 4: 60%


      bucket 2: 40%


      bucket 3: 60%


      bucket 6: 40%


      bucket 8: 40%


      bucket 1: 70%


      bucket 4: 70%


      bucket 3: 70%


      bucket 5: 50%


      bucket 2: 50%


      bucket 7: 30%


      bucket 6: 50%


      bucket 1: 80%


      bucket 8: 50%


      bucket 4: 80%


      bucket 3: 80%


      bucket 5: 60%


      bucket 1: 90%


      bucket 2: 60%


      bucket 6: 60%


      bucket 8: 60%


      bucket 4: 90%


      bucket 3: 90%


      bucket 7: 40%


      bucket 1: 100%


      Sorting block of length 854 for bucket 1


      (Using difference cover)


      bucket 5: 70%


      bucket 2: 70%


      bucket 4: 100%


      Sorting block of length 602 for bucket 4


      (Using difference cover)


      bucket 3: 100%


      bucket 8: 70%


      Sorting block of length 814 for bucket 3


      (Using difference cover)


      bucket 6: 70%


      bucket 5: 80%


      bucket 2: 80%


      bucket 7: 50%


      bucket 8: 80%


      bucket 6: 80%


      bucket 2: 90%


      bucket 5: 90%


      bucket 8: 90%


      bucket 2: 100%


      Sorting block of length 671 for bucket 2


      (Using difference cover)


      bucket 6: 90%


      Sorting block time: 00:00:00


      bucket 7: 60%


      bucket 5: 100%


    Returning block of 815 for bucket 3


      bucket 8: 100%


      Sorting block of length 754 for bucket 5


      (Using difference cover)


      Sorting block of length 902 for bucket 8


      (Using difference cover)


      Sorting block time: 00:00:00


    Returning block of 603 for bucket 4


      bucket 6: 100%


      Sorting block of length 214 for bucket 6


      (Using difference cover)


      bucket 7: 70%


    Getting block 9 of 60


      Reserving size (936) for bucket 9


      Calculating Z arrays for bucket 9


      Entering block accumulator loop for bucket 9:


      Sorting block time: 00:00:00


    Returning block of 855 for bucket 1


      Sorting block time: 00:00:00


    Returning block of 755 for bucket 5


    Getting block 10 of 60


      bucket 9: 10%


      Reserving size (936) for bucket 10


      Calculating Z arrays for bucket 10


      bucket 7: 80%


      Entering block accumulator loop for bucket 10:


      bucket 9: 20%


      Sorting block time: 00:00:00


    Returning block of 903 for bucket 8


      bucket 10: 10%


      Sorting block time: 00:00:00


      bucket 9: 30%


    Returning block of 672 for bucket 2


    Getting block 11 of 60


      Reserving size (936) for bucket 11


    Getting block 12 of 60


      Calculating Z arrays for bucket 11


      bucket 7: 90%


      Reserving size (936) for bucket 12


      bucket 10: 20%


      Sorting block time: 00:00:00


      bucket 9: 40%


      Entering block accumulator loop for bucket 11:


      Calculating Z arrays for bucket 12


    Getting block 13 of 60


    Returning block of 215 for bucket 6


      Reserving size (936) for bucket 13


      Entering block accumulator loop for bucket 12:


      Calculating Z arrays for bucket 13


    Getting block 14 of 60


      Reserving size (936) for bucket 14


      Calculating Z arrays for bucket 14


      Entering block accumulator loop for bucket 13:


      bucket 9: 50%


      Entering block accumulator loop for bucket 14:


      bucket 10: 30%


      bucket 7: 100%


      bucket 11: 10%


      bucket 12: 10%


      Sorting block of length 914 for bucket 7


      (Using difference cover)


      bucket 9: 60%


      bucket 13: 10%


    Getting block 15 of 60


      bucket 10: 40%


      Reserving size (936) for bucket 15


      bucket 14: 10%


      Calculating Z arrays for bucket 15


      bucket 11: 20%


      bucket 12: 20%


      bucket 9: 70%


      Entering block accumulator loop for bucket 15:


      bucket 13: 20%


      bucket 10: 50%


      bucket 14: 20%


      bucket 9: 80%


      bucket 11: 30%


      bucket 12: 30%


      bucket 13: 30%


      bucket 10: 60%


      bucket 15: 10%


      bucket 14: 30%


      bucket 9: 90%


      bucket 11: 40%


      bucket 12: 40%


      bucket 13: 40%


      bucket 10: 70%


      bucket 9: 100%


      Sorting block of length 908 for bucket 9


      (Using difference cover)


      bucket 14: 40%


      bucket 11: 50%


      bucket 12: 50%


      bucket 13: 50%


      bucket 10: 80%


      Sorting block time: 00:00:00


      bucket 15: 20%


      bucket 14: 50%


    Returning block of 915 for bucket 7


      bucket 12: 60%


      Sorting block time: 00:00:00


      bucket 10: 90%


      bucket 11: 60%


      bucket 13: 60%


      bucket 14: 60%


    Returning block of 909 for bucket 9


      bucket 15: 30%


      bucket 10: 100%


      Sorting block of length 806 for bucket 10


      (Using difference cover)


      bucket 12: 70%


      bucket 13: 70%


      bucket 11: 70%


      bucket 14: 70%


      bucket 15: 40%


    Getting block 16 of 60


      Reserving size (936) for bucket 16


      Calculating Z arrays for bucket 16


      bucket 12: 80%


      Entering block accumulator loop for bucket 16:


      bucket 13: 80%


      bucket 14: 80%


    Getting block 17 of 60


      bucket 11: 80%


      Sorting block time: 00:00:00


      Reserving size (936) for bucket 17


    Returning block of 807 for bucket 10


      Calculating Z arrays for bucket 17


      bucket 15: 50%


      bucket 12: 90%


      bucket 16: 10%


      bucket 13: 90%


      Entering block accumulator loop for bucket 17:


      bucket 14: 90%


      bucket 11: 90%


    Getting block 18 of 60


      bucket 12: 100%


      bucket 16: 20%


      bucket 13: 100%


      Reserving size (936) for bucket 18


      Sorting block of length 684 for bucket 12


      (Using difference cover)


      bucket 14: 100%


      bucket 15: 60%


      Sorting block of length 290 for bucket 13


      (Using difference cover)


      Calculating Z arrays for bucket 18


      Sorting block of length 759 for bucket 14


      (Using difference cover)


      bucket 11: 100%


      bucket 17: 10%


      Sorting block of length 538 for bucket 11


      (Using difference cover)


      Entering block accumulator loop for bucket 18:


      bucket 16: 30%


      bucket 15: 70%


      bucket 16: 40%


      bucket 18: 10%


      bucket 17: 20%


      Sorting block time: 00:00:00


    Returning block of 685 for bucket 12


      Sorting block time: 00:00:00


    Returning block of 291 for bucket 13


      bucket 16: 50%


      bucket 18: 20%


      Sorting block time: 00:00:00


    Returning block of 539 for bucket 11


      bucket 15: 80%


      bucket 17: 30%


    Getting block 19 of 60


      Sorting block time: 00:00:00


    Getting block 20 of 60


      Reserving size (936) for bucket 19


    Returning block of 760 for bucket 14


      Reserving size (936) for bucket 20


      bucket 16: 60%


      Calculating Z arrays for bucket 19


      bucket 18: 30%


      Calculating Z arrays for bucket 20


    Getting block 21 of 60


      Entering block accumulator loop for bucket 19:


      Reserving size (936) for bucket 21


      Calculating Z arrays for bucket 21


      Entering block accumulator loop for bucket 20:


    Getting block 22 of 60


      Reserving size (936) for bucket 22


      Entering block accumulator loop for bucket 21:


      Calculating Z arrays for bucket 22


      Entering block accumulator loop for bucket 22:


      bucket 18: 40%


      bucket 16: 70%


      bucket 15: 90%


      bucket 17: 40%


      bucket 19: 10%


      bucket 20: 10%


      bucket 18: 50%


      bucket 21: 10%


      bucket 22: 10%


      bucket 16: 80%


      bucket 19: 20%


      bucket 18: 60%


      bucket 17: 50%


      bucket 15: 100%


      bucket 20: 20%


      Sorting block of length 715 for bucket 15


      (Using difference cover)


      bucket 16: 90%


      bucket 21: 20%


      bucket 22: 20%


      bucket 18: 70%


      bucket 19: 30%


      bucket 17: 60%


      bucket 18: 80%


      bucket 16: 100%


      bucket 20: 30%


      Sorting block of length 865 for bucket 16


      (Using difference cover)


      bucket 22: 30%


      bucket 21: 30%


      bucket 19: 40%


      bucket 18: 90%


      Sorting block time: 00:00:00


    Returning block of 716 for bucket 15


      bucket 20: 40%


      bucket 17: 70%


      Sorting block time: 00:00:00


    Returning block of 866 for bucket 16


      bucket 22: 40%


      bucket 21: 40%


      bucket 18: 100%


      Sorting block of length 604 for bucket 18


      (Using difference cover)


      bucket 19: 50%


      bucket 20: 50%


    Getting block 23 of 60


      bucket 22: 50%


      Reserving size (936) for bucket 23


      bucket 21: 50%


      Calculating Z arrays for bucket 23


      bucket 17: 80%


    Getting block 24 of 60


      bucket 19: 60%


      Reserving size (936) for bucket 24


      Sorting block time: 00:00:00


      Entering block accumulator loop for bucket 23:


      bucket 20: 60%


      Calculating Z arrays for bucket 24


    Returning block of 605 for bucket 18


      Entering block accumulator loop for bucket 24:


      bucket 22: 60%


      bucket 21: 60%


    Getting block 25 of 60


      Reserving size (936) for bucket 25


      Calculating Z arrays for bucket 25


      Entering block accumulator loop for bucket 25:


      bucket 19: 70%


      bucket 24: 10%


      bucket 20: 70%


      bucket 17: 90%


      bucket 22: 70%


      bucket 21: 70%


      bucket 25: 10%


      bucket 23: 10%


      bucket 24: 20%


      bucket 19: 80%


      bucket 20: 80%


      bucket 22: 80%


      bucket 25: 20%


      bucket 21: 80%


      bucket 17: 100%


      Sorting block of length 817 for bucket 17


      (Using difference cover)


      bucket 24: 30%


      bucket 19: 90%


      bucket 20: 90%


      bucket 22: 90%


      bucket 23: 20%


      bucket 25: 30%


      bucket 21: 90%


      bucket 24: 40%


      bucket 20: 100%


      bucket 19: 100%


      bucket 22: 100%


      Sorting block of length 408 for bucket 20


      Sorting block time: 00:00:00


      (Using difference cover)


      bucket 24: 50%


      bucket 21: 100%


      bucket 25: 40%


      Sorting block of length 736 for bucket 19


      (Using difference cover)


      Sorting block of length 839 for bucket 22


      (Using difference cover)


    Returning block of 818 for bucket 17


      Sorting block of length 681 for bucket 21


      (Using difference cover)


      bucket 23: 30%


      bucket 24: 60%


      bucket 25: 50%


      Sorting block time: 00:00:00


    Returning block of 409 for bucket 20


    Getting block 26 of 60


      Sorting block time: 00:00:00


      Reserving size (936) for bucket 26


    Returning block of 840 for bucket 22


      Calculating Z arrays for bucket 26


      Sorting block time: 00:00:00


      bucket 24: 70%


    Getting block 27 of 60


    Returning block of 682 for bucket 21


      Sorting block time: 00:00:00


      bucket 23: 40%


      bucket 25: 60%


      Reserving size (936) for bucket 27


      Entering block accumulator loop for bucket 26:


    Returning block of 737 for bucket 19


      bucket 24: 80%


      Calculating Z arrays for bucket 27


    Getting block 28 of 60


    Getting block 29 of 60


      Reserving size (936) for bucket 28


      Reserving size (936) for bucket 29


      Calculating Z arrays for bucket 28


      Entering block accumulator loop for bucket 27:


      Calculating Z arrays for bucket 29


      Entering block accumulator loop for bucket 28:


      Entering block accumulator loop for bucket 29:


      bucket 25: 70%


    Getting block 30 of 60


      Reserving size (936) for bucket 30


      Calculating Z arrays for bucket 30


      Entering block accumulator loop for bucket 30:


      bucket 24: 90%


      bucket 23: 50%


      bucket 25: 80%


      bucket 27: 10%


      bucket 26: 10%


      bucket 28: 10%


      bucket 29: 10%


      bucket 30: 10%


      bucket 25: 90%


      bucket 24: 100%


      Sorting block of length 461 for bucket 24


      (Using difference cover)


      bucket 28: 20%


      bucket 27: 20%


      bucket 25: 100%


      Sorting block of length 879 for bucket 25


      (Using difference cover)


      bucket 26: 20%


      bucket 23: 60%


      bucket 29: 20%


      bucket 30: 20%


      bucket 27: 30%


      Sorting block time: 00:00:00


      bucket 28: 30%


      bucket 26: 30%


      Sorting block time: 00:00:00


    Returning block of 462 for bucket 24


    Returning block of 880 for bucket 25


      bucket 29: 30%


    Getting block 31 of 60


      Reserving size (936) for bucket 31


      bucket 23: 70%


      Calculating Z arrays for bucket 31


      bucket 30: 30%


      Entering block accumulator loop for bucket 31:


    Getting block 32 of 60


      bucket 27: 40%


      Reserving size (936) for bucket 32


      bucket 28: 40%


      Calculating Z arrays for bucket 32


      Entering block accumulator loop for bucket 32:


      bucket 26: 40%


      bucket 29: 40%


      bucket 31: 10%


      bucket 30: 40%


      bucket 23: 80%


      bucket 27: 50%


      bucket 32: 10%


      bucket 28: 50%


      bucket 31: 20%


      bucket 26: 50%


      bucket 30: 50%


      bucket 29: 50%


      bucket 32: 20%


      bucket 27: 60%


      bucket 28: 60%


      bucket 32: 30%


      bucket 31: 30%


      bucket 23: 90%


      bucket 30: 60%


      bucket 26: 60%


      bucket 27: 70%


      bucket 29: 60%


      bucket 32: 40%


      bucket 28: 70%


      bucket 31: 40%


      bucket 30: 70%


      bucket 32: 50%


      bucket 27: 80%


      bucket 23: 100%


      bucket 31: 50%


      Sorting block of length 601 for bucket 23


      bucket 28: 80%


      (Using difference cover)


      bucket 29: 70%


      bucket 26: 70%


      bucket 30: 80%


      bucket 32: 60%


      bucket 27: 90%


      bucket 31: 60%


      bucket 28: 90%


      bucket 32: 70%


      bucket 29: 80%


      bucket 30: 90%


      bucket 31: 70%


      bucket 26: 80%


      bucket 27: 100%


      Sorting block of length 827 for bucket 27


      (Using difference cover)


      bucket 32: 80%


      bucket 28: 100%


      Sorting block of length 862 for bucket 28


      (Using difference cover)


      bucket 31: 80%


      bucket 30: 100%


      bucket 29: 90%


      Sorting block of length 449 for bucket 30


      (Using difference cover)


      bucket 32: 90%


      bucket 31: 90%


      bucket 26: 90%


      Sorting block time: 00:00:00


    Returning block of 602 for bucket 23


      bucket 32: 100%


      bucket 29: 100%


      Sorting block of length 607 for bucket 32


      Sorting block time: 00:00:00


      (Using difference cover)


      Sorting block of length 616 for bucket 29


      (Using difference cover)


    Returning block of 828 for bucket 27


      bucket 31: 100%


      Sorting block of length 720 for bucket 31


      (Using difference cover)


      bucket 26: 100%


      Sorting block of length 732 for bucket 26


      Sorting block time: 00:00:00


    Returning block of 450 for bucket 30


    Getting block 33 of 60


      (Using difference cover)


    Getting block 34 of 60


      Reserving size (936) for bucket 33


      Reserving size (936) for bucket 34


      Sorting block time: 00:00:00


      Calculating Z arrays for bucket 33


      Calculating Z arrays for bucket 34


    Returning block of 863 for bucket 28


      Entering block accumulator loop for bucket 33:


    Getting block 35 of 60


      Reserving size (936) for bucket 35


      Sorting block time: 00:00:00


      Calculating Z arrays for bucket 35


      Entering block accumulator loop for bucket 34:


    Returning block of 608 for bucket 32


      Entering block accumulator loop for bucket 35:


      bucket 33: 10%


      Sorting block time: 00:00:00


    Getting block 36 of 60


    Returning block of 617 for bucket 29


      Reserving size (936) for bucket 36


      Calculating Z arrays for bucket 36


      Sorting block time: 00:00:00


    Returning block of 721 for bucket 31


    Getting block 37 of 60


      Entering block accumulator loop for bucket 36:


      Reserving size (936) for bucket 37


      bucket 35: 10%


      Calculating Z arrays for bucket 37


      Entering block accumulator loop for bucket 37:


      bucket 33: 20%


    Getting block 38 of 60


      Reserving size (936) for bucket 38


      Calculating Z arrays for bucket 38


      bucket 34: 10%


      Entering block accumulator loop for bucket 38:


    Getting block 39 of 60


      bucket 36: 10%


      Reserving size (936) for bucket 39


      bucket 35: 20%


      Calculating Z arrays for bucket 39


      bucket 37: 10%


      Sorting block time: 00:00:00


    Returning block of 733 for bucket 26


      bucket 33: 30%


      Entering block accumulator loop for bucket 39:


      bucket 38: 10%


      bucket 36: 20%


      bucket 37: 20%


      bucket 35: 30%


      bucket 34: 20%


      bucket 39: 10%


      bucket 33: 40%


    Getting block 40 of 60


      Reserving size (936) for bucket 40


      Calculating Z arrays for bucket 40


      Entering block accumulator loop for bucket 40:


      bucket 38: 20%


      bucket 37: 30%


      bucket 36: 30%


      bucket 39: 20%


      bucket 35: 40%


      bucket 33: 50%


      bucket 34: 30%


      bucket 37: 40%


      bucket 38: 30%


      bucket 36: 40%


      bucket 39: 30%


      bucket 35: 50%


      bucket 33: 60%


      bucket 40: 10%


      bucket 37: 50%


      bucket 38: 40%


      bucket 36: 50%


      bucket 39: 40%


      bucket 34: 40%


      bucket 35: 60%


      bucket 33: 70%


      bucket 40: 20%


      bucket 37: 60%


      bucket 39: 50%


      bucket 36: 60%


      bucket 38: 50%


      bucket 35: 70%


      bucket 33: 80%


      bucket 34: 50%


      bucket 39: 60%


      bucket 37: 70%


      bucket 40: 30%


      bucket 36: 70%


      bucket 38: 60%


      bucket 35: 80%


      bucket 33: 90%


      bucket 39: 70%


      bucket 37: 80%


      bucket 36: 80%


      bucket 34: 60%


      bucket 40: 40%


      bucket 35: 90%


      bucket 38: 70%


      bucket 33: 100%


      bucket 39: 80%


      Sorting block of length 608 for bucket 33


      (Using difference cover)


      bucket 37: 90%


      bucket 36: 90%


      bucket 40: 50%


      bucket 39: 90%


      bucket 38: 80%


      bucket 35: 100%


      Sorting block time: 00:00:00


      bucket 34: 70%


      Sorting block of length 346 for bucket 35


      (Using difference cover)


    Returning block of 609 for bucket 33


      bucket 40: 60%


      bucket 36: 100%


      bucket 37: 100%


      bucket 39: 100%


      Sorting block of length 446 for bucket 36


      (Using difference cover)


      bucket 38: 90%


      Sorting block of length 926 for bucket 37


      (Using difference cover)


      Sorting block of length 988 for bucket 39


      (Using difference cover)


    Getting block 41 of 60


      Reserving size (936) for bucket 41


      Calculating Z arrays for bucket 41


      Entering block accumulator loop for bucket 41:


      Sorting block time: 00:00:00


    Returning block of 347 for bucket 35


      bucket 40: 70%


      bucket 34: 80%


      bucket 38: 100%


      Sorting block of length 774 for bucket 38


      (Using difference cover)


      Sorting block time: 00:00:00


    Returning block of 447 for bucket 36


    Getting block 42 of 60


      bucket 41: 10%


      Reserving size (936) for bucket 42


      Calculating Z arrays for bucket 42


      Entering block accumulator loop for bucket 42:


      bucket 40: 80%


    Getting block 43 of 60


      Reserving size (936) for bucket 43


      Calculating Z arrays for bucket 43


      Entering block accumulator loop for bucket 43:


      Sorting block time: 00:00:00


    Returning block of 927 for bucket 37


      bucket 41: 20%


      bucket 34: 90%


      bucket 40: 90%


      bucket 42: 10%


      Sorting block time: 00:00:00


    Returning block of 775 for bucket 38


      bucket 43: 10%


    Getting block 44 of 60


      bucket 40: 100%


      bucket 41: 30%


      Reserving size (936) for bucket 44


      Sorting block of length 146 for bucket 40


      (Using difference cover)


      Calculating Z arrays for bucket 44


    Getting block 45 of 60


      Reserving size (936) for bucket 45


      Calculating Z arrays for bucket 45


      bucket 42: 20%


      Sorting block time: 00:00:00


      Entering block accumulator loop for bucket 44:


    Returning block of 989 for bucket 39


      Entering block accumulator loop for bucket 45:


      bucket 34: 100%


      bucket 43: 20%


      Sorting block of length 631 for bucket 34


      (Using difference cover)


      bucket 41: 40%


      Sorting block time: 00:00:00


    Returning block of 147 for bucket 40


      bucket 42: 30%


    Getting block 46 of 60


      Reserving size (936) for bucket 46


      Calculating Z arrays for bucket 46


    Getting block 47 of 60


      bucket 44: 10%


      Reserving size (936) for bucket 47


      Entering block accumulator loop for bucket 46:


      Calculating Z arrays for bucket 47


      bucket 43: 30%


      Entering block accumulator loop for bucket 47:


      bucket 41: 50%


      bucket 45: 10%


      Sorting block time: 00:00:00


    Returning block of 632 for bucket 34


      bucket 46: 10%


      bucket 42: 40%


      bucket 44: 20%


      bucket 43: 40%


      bucket 41: 60%


      bucket 45: 20%


      bucket 46: 20%


      bucket 44: 30%


      bucket 47: 10%


      bucket 42: 50%


    Getting block 48 of 60


      Reserving size (936) for bucket 48


      Calculating Z arrays for bucket 48


      bucket 43: 50%


      bucket 46: 30%


      Entering block accumulator loop for bucket 48:


      bucket 44: 40%


      bucket 41: 70%


      bucket 47: 20%


      bucket 45: 30%


      bucket 46: 40%


      bucket 42: 60%


      bucket 43: 60%


      bucket 44: 50%


      bucket 47: 30%


      bucket 41: 80%


      bucket 46: 50%


      bucket 48: 10%


      bucket 44: 60%


      bucket 42: 70%


      bucket 47: 40%


      bucket 45: 40%


      bucket 43: 70%


      bucket 41: 90%


      bucket 46: 60%


      bucket 47: 50%


      bucket 44: 70%


      bucket 42: 80%


      bucket 43: 80%


      bucket 48: 20%


      bucket 45: 50%


      bucket 46: 70%


      bucket 41: 100%


      bucket 47: 60%


      bucket 44: 80%


      Sorting block of length 554 for bucket 41


      (Using difference cover)


      bucket 42: 90%


      bucket 43: 90%


      bucket 46: 80%


      bucket 47: 70%


      bucket 45: 60%


      bucket 44: 90%


      bucket 48: 30%


      bucket 46: 90%


      bucket 42: 100%


      bucket 47: 80%


      Sorting block of length 478 for bucket 42


      (Using difference cover)


      bucket 43: 100%


      bucket 44: 100%


      Sorting block time: 00:00:00


      Sorting block of length 734 for bucket 43


      (Using difference cover)


      Sorting block of length 628 for bucket 44


      (Using difference cover)


    Returning block of 555 for bucket 41


      bucket 45: 70%


      bucket 46: 100%


      bucket 47: 90%


      bucket 48: 40%


      Sorting block of length 683 for bucket 46


      (Using difference cover)


    Getting block 49 of 60


      Reserving size (936) for bucket 49


      Calculating Z arrays for bucket 49


      bucket 45: 80%


      Entering block accumulator loop for bucket 49:


      Sorting block time: 00:00:00


    Returning block of 479 for bucket 42


      bucket 47: 100%


      Sorting block of length 434 for bucket 47


      (Using difference cover)


      Sorting block time: 00:00:00


    Returning block of 735 for bucket 43


      bucket 48: 50%


    Getting block 50 of 60


      Reserving size (936) for bucket 50


      Calculating Z arrays for bucket 50


      bucket 45: 90%


      Sorting block time: 00:00:00


      bucket 49: 10%


    Returning block of 629 for bucket 44


      Entering block accumulator loop for bucket 50:


      Sorting block time: 00:00:00


    Returning block of 684 for bucket 46


    Getting block 51 of 60


      Reserving size (936) for bucket 51


      Calculating Z arrays for bucket 51


      Entering block accumulator loop for bucket 51:


    Getting block 52 of 60


      bucket 48: 60%


      bucket 49: 20%


      bucket 50: 10%


      bucket 45: 100%


      Sorting block time: 00:00:00


      Reserving size (936) for bucket 52


      Sorting block of length 883 for bucket 45


      (Using difference cover)


    Returning block of 435 for bucket 47


      Calculating Z arrays for bucket 52


    Getting block 53 of 60


      Reserving size (936) for bucket 53


      Calculating Z arrays for bucket 53


      bucket 51: 10%


      Entering block accumulator loop for bucket 52:


      Entering block accumulator loop for bucket 53:


      bucket 49: 30%


    Getting block 54 of 60


      bucket 50: 20%


      Reserving size (936) for bucket 54


      bucket 48: 70%


      Calculating Z arrays for bucket 54


      bucket 52: 10%


      Entering block accumulator loop for bucket 54:


      bucket 51: 20%


      bucket 53: 10%


      bucket 49: 40%


      Sorting block time: 00:00:00


    Returning block of 884 for bucket 45


      bucket 50: 30%


      bucket 52: 20%


      bucket 54: 10%


      bucket 53: 20%


      bucket 48: 80%


      bucket 51: 30%


      bucket 49: 50%


    Getting block 55 of 60


      Reserving size (936) for bucket 55


      Calculating Z arrays for bucket 55


      bucket 52: 30%


      bucket 50: 40%


      Entering block accumulator loop for bucket 55:


      bucket 53: 30%


      bucket 54: 20%


      bucket 49: 60%


      bucket 51: 40%


      bucket 48: 90%


      bucket 53: 40%


      bucket 52: 40%


      bucket 54: 30%


      bucket 50: 50%


      bucket 55: 10%


      bucket 53: 50%


      bucket 49: 70%


      bucket 51: 50%


      bucket 52: 50%


      bucket 48: 100%


      Sorting block of length 643 for bucket 48


      (Using difference cover)


      bucket 54: 40%


      bucket 53: 60%


      bucket 50: 60%


      bucket 55: 20%


      bucket 49: 80%


      bucket 51: 60%


      bucket 52: 60%


      bucket 53: 70%


      bucket 54: 50%


      bucket 50: 70%


      Sorting block time: 00:00:00


    Returning block of 644 for bucket 48


      bucket 49: 90%


      bucket 52: 70%


      bucket 55: 30%


      bucket 53: 80%


      bucket 51: 70%


      bucket 54: 60%


      bucket 52: 80%


      bucket 50: 80%


      bucket 53: 90%


      bucket 49: 100%


    Getting block 56 of 60


      Sorting block of length 757 for bucket 49


      (Using difference cover)


      Reserving size (936) for bucket 56


      Calculating Z arrays for bucket 56


      bucket 55: 40%


      bucket 51: 80%


      Entering block accumulator loop for bucket 56:


      bucket 54: 70%


      bucket 52: 90%


      bucket 53: 100%


      Sorting block of length 614 for bucket 53


      (Using difference cover)


      bucket 50: 90%


      bucket 55: 50%


      bucket 52: 100%


      Sorting block of length 613 for bucket 52


      Sorting block time: 00:00:00


      (Using difference cover)


      bucket 51: 90%


    Returning block of 758 for bucket 49


      bucket 54: 80%


      Sorting block time: 00:00:00


    Returning block of 615 for bucket 53


      bucket 56: 10%


      bucket 55: 60%


      bucket 50: 100%


    Getting block 57 of 60


      Sorting block time: 00:00:00


      Sorting block of length 302 for bucket 50


      (Using difference cover)


      Reserving size (936) for bucket 57


    Returning block of 614 for bucket 52


      Calculating Z arrays for bucket 57


      bucket 51: 100%


      Sorting block of length 653 for bucket 51


      (Using difference cover)


      bucket 54: 90%


      Entering block accumulator loop for bucket 57:


    Getting block 58 of 60


      Reserving size (936) for bucket 58


      Calculating Z arrays for bucket 58


      Entering block accumulator loop for bucket 58:


    Getting block 59 of 60


      bucket 55: 70%


      Reserving size (936) for bucket 59


      Sorting block time: 00:00:00


      Sorting block time: 00:00:00


    Returning block of 654 for bucket 51


    Returning block of 303 for bucket 50


      Calculating Z arrays for bucket 59


      bucket 54: 100%


      bucket 56: 20%


      bucket 57: 10%


      Sorting block of length 507 for bucket 54


      (Using difference cover)


      Entering block accumulator loop for bucket 59:


      bucket 58: 10%


      bucket 55: 80%


    Getting block 60 of 60


      Reserving size (936) for bucket 60


      Calculating Z arrays for bucket 60


      Entering block accumulator loop for bucket 60:


      bucket 60: 10%


      bucket 58: 20%


      bucket 57: 20%


      Sorting block time: 00:00:00


    Returning block of 508 for bucket 54


      bucket 59: 10%


      bucket 56: 30%


      bucket 60: 20%


      bucket 55: 90%


      bucket 58: 30%


      bucket 57: 30%


      bucket 60: 30%


      bucket 58: 40%


      bucket 59: 20%


      bucket 60: 40%


      bucket 57: 40%


      bucket 55: 100%


      Sorting block of length 812 for bucket 55


      (Using difference cover)


      bucket 58: 50%


      bucket 59: 30%


      bucket 60: 50%


      bucket 56: 40%


      bucket 58: 60%


      bucket 57: 50%


      bucket 60: 60%


      bucket 59: 40%


      Sorting block time: 00:00:00


    Returning block of 813 for bucket 55


      bucket 58: 70%


      bucket 60: 70%


      bucket 59: 50%


      bucket 57: 60%


      bucket 56: 50%


      bucket 60: 80%


      bucket 58: 80%


      bucket 59: 60%


      bucket 57: 70%


      bucket 60: 90%


      bucket 59: 70%


      bucket 58: 90%


      bucket 57: 80%


      bucket 60: 100%


      Sorting block of length 754 for bucket 60


      (Using difference cover)


      bucket 56: 60%


      bucket 58: 100%


      Sorting block of length 603 for bucket 58


      (Using difference cover)


      bucket 59: 80%


      bucket 57: 90%


      Sorting block time: 00:00:00


    Returning block of 755 for bucket 60


      bucket 59: 90%


      bucket 57: 100%


      Sorting block of length 633 for bucket 57


      (Using difference cover)


      Sorting block time: 00:00:00


      bucket 56: 70%


    Returning block of 604 for bucket 58


      bucket 59: 100%


      Sorting block of length 878 for bucket 59


      (Using difference cover)


      Sorting block time: 00:00:00


      bucket 56: 80%


    Returning block of 634 for bucket 57


      Sorting block time: 00:00:00


    Returning block of 879 for bucket 59


      bucket 56: 90%


      bucket 56: 100%


      Sorting block of length 721 for bucket 56


      (Using difference cover)


      Sorting block time: 00:00:00


    Returning block of 722 for bucket 56


    Exited Ebwt loop


    fchr[A]: 0


    fchr[C]: 10842


    fchr[G]: 19880


    fchr[T]: 30171


    fchr[$]: 39937


    Exiting Ebwt::buildToDisk()


    Returning from initFromVector


    Wrote 4207850 bytes to primary EBWT file: /tmp/refgenie_demo_2hzk45es/data/kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB/bowtie2_index/2.3.5/kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB.rev.1.bt2


    Wrote 9992 bytes to secondary EBWT file: /tmp/refgenie_demo_2hzk45es/data/kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB/bowtie2_index/2.3.5/kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB.rev.2.bt2


    Re-opening _in1 and _in2 as input streams


    Returning from Ebwt constructor


    Headers:


        len: 39937


        bwtLen: 39938


        sz: 9985


        bwtSz: 9985


        lineRate: 6


        offRate: 4


        offMask: 0xfffffff0


        ftabChars: 10


        eftabLen: 20


        eftabSz: 80


        ftabLen: 1048577


        ftabSz: 4194308


        offsLen: 2497


        offsSz: 9988


        lineSz: 64


        sideSz: 64


        sideBwtSz: 48


        sideBwtLen: 192


        numSides: 209


        numLines: 209


        ebwtTotLen: 13376


        ebwtTotSz: 13376


        color: 0


        reverse: 1


    Total time for backward call to driver() for mirror index: 00:00:00


    </pre>


    Command completed. Elapsed time: 0:00:01. Running peak memory: 0.037GB.  
      PID: 1045591;	Command: bowtie2-build;	Return code: 0;	Memory used: 0.037GB
    



<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Asset <span style="color: #008000; text-decoration-color: #008000">'t7/bowtie2_index:2.3.5'</span> build succeeded                                              <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/builder.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">builder.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/builder.py#666" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">666</span></a>
</pre>



    
    ### Pipeline completed. Epilogue


    *        Elapsed time (this run):  0:00:01


    *  Total elapsed time (all runs):  0:00:00


    *         Peak memory (this run):  0.0368 GB


    *        Pipeline completed time: 2026-08-20 21:55:45



<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"></pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Added: <span style="color: #008000; text-decoration-color: #008000">'t7/bowtie2_index:2.3.5'</span>                                                             <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/content.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">content.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/content.py#381" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">381</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Set default asset: <span style="color: #008000; text-decoration-color: #008000">'kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB/bowtie2_index:2.3.5'</span>                      <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/seek.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">seek.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/seek.py#612" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">612</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Added asset: <span style="color: #008000; text-decoration-color: #008000">'kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB/bowtie2_index:2.3.5'</span>                         <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/builder.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">builder.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/builder.py#705" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">705</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Created alias directories: <span style="color: #800080; text-decoration-color: #800080">/tmp/refgenie_demo_2hzk45es/alias/t7/bowtie2_index/</span><span style="color: #ff00ff; text-decoration-color: #ff00ff">2.3.5</span>        <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/symlinks.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">symlinks.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/symlinks.py#223" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">223</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Updated parents of <span style="color: #008000; text-decoration-color: #008000">'kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB/bowtie2_index:2.3.5'</span>                 <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/relations.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">relations.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/relations.py#166" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">166</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Staging asset: kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB/bowtie2_index:<span style="color: #008080; text-decoration-color: #008080; font-weight: bold">2.3</span>.<span style="color: #008080; text-decoration-color: #008080; font-weight: bold">5</span>                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/stage.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">stage.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/stage.py#67" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">67</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reconciled file-serving stage dir:                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/stage.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">stage.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/stage.py#335" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">335</span></a>
         <span style="color: #800080; text-decoration-color: #800080">/tmp/refgenie_archive_demo_3ps6xfvv/kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB/bowtie2_index/</span><span style="color: #ff00ff; text-decoration-color: #ff00ff">2.3.5</span> <span style="font-weight: bold">(</span><span style="color: #008080; text-decoration-color: #008080; font-weight: bold">6</span>   <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">            </span>
         <span style="color: #800080; text-decoration-color: #800080; font-weight: bold">link</span><span style="font-weight: bold">(</span>s<span style="font-weight: bold">)</span>, <span style="color: #008080; text-decoration-color: #008080; font-weight: bold">6</span> written or removed; excluded <span style="color: #008080; text-decoration-color: #008080; font-weight: bold">0</span> colocation <span style="color: #800080; text-decoration-color: #800080; font-weight: bold">symlink</span><span style="font-weight: bold">(</span>s<span style="font-weight: bold">))</span>                              <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">            </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Created file-serving StagedAsset for kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB/bowtie2_index:<span style="color: #008080; text-decoration-color: #008080; font-weight: bold">2.3</span>.<span style="color: #008080; text-decoration-color: #008080; font-weight: bold">5</span>     <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/stage.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">stage.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/stage.py#217" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">217</span></a>
</pre>






    Asset(name='2.3.5', description='Genome index for bowtie2, produced with bowtie2-build', size=8445686, colocate=None, updated_at=datetime.datetime(2026, 8, 21, 1, 55, 46, 9717), path='data/kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB/bowtie2_index/2b4c1900ec5a2b9fce39e45a71e6bed38ca234d394b8af28c2a3b0e09548d51e', digest='2b4c1900ec5a2b9fce39e45a71e6bed38ca234d394b8af28c2a3b0e09548d51e', recipe_id=2, serving_modes_override=None, asset_group_id=3, created_at=datetime.datetime(2026, 8, 21, 1, 55, 46, 9733))



Let's list the assets for the genome `t7` to verify that the `bowtie2_index` asset has been built successfully.



```python
refgenie.asset.table(genome_digests=[refgenie.alias.resolve(GenomeAlias("t7"))])[0]
```




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="font-style: italic">                     Refgenie assets. Source: local                     </span>
┏━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━┳━━━━━━━━━┓
┃<span style="font-weight: bold"> Aliases </span>┃<span style="font-weight: bold"> Genome digest                    </span>┃<span style="font-weight: bold"> Asset group   </span>┃<span style="font-weight: bold"> Asset   </span>┃
┡━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━╇━━━━━━━━━┩
│ t7      │ kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB │ fasta         │ default │
│ t7      │ kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB │ bowtie2_index │ 2.3.5   │
└─────────┴──────────────────────────────────┴───────────────┴─────────┘
</pre>




One of the assets was also archived (a neccessary step to serve the assets via the refgenie server). Let's list the archived assets.


```python
print(refgenie.stage.table())
```


<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="font-style: italic">                                                   Staged Assets                                                   </span>
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━┳━━━━━━┓
┃<span style="font-weight: bold"> Asset Digest                                    </span>┃<span style="font-weight: bold"> Asset name                                      </span>┃<span style="font-weight: bold"> Mode </span>┃<span style="font-weight: bold"> Size </span>┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━╇━━━━━━┩
│ 2b4c1900ec5a2b9fce39e45a71e6bed38ca234d394b8af… │ kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB/bowtie2_index… │ file │ -    │
└─────────────────────────────────────────────────┴─────────────────────────────────────────────────┴──────┴──────┘
</pre>



Asset `bowtie2_index` has been built successfully for the `t7` genome, and automatically tagged with `2.3.5`, indicating the version of Bowtie2 software used (this behavior is encoded in the recipe).

## Interact with aliases

Let's list the aliases:


```python
refgenie.alias.table()
```




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="font-style: italic">                Genome aliases                </span>
┏━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃<span style="font-weight: bold"> Aliases </span>┃<span style="font-weight: bold"> Genome digest                    </span>┃
┡━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ rCRSd   │ jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP │
│ t7      │ kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB │
└─────────┴──────────────────────────────────┘
</pre>




Let's assign another alias to the same genome digest, this way we can refer to the same genome in multiple ways.


```python
t7_alias = refgenie.set_genome_alias(
    alias_name="Bacteriophage-T7",
    genome_digest="kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB",
    genome_description="My favorite genome",
)
fav_alias = refgenie.set_genome_alias(
    alias_name="myFavGenome",
    genome_digest="kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB",
    genome_description="My favorite genome",
)
```


<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Added alias: Bacteriophage-T7                                                                 <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/alias.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">alias.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/alias.py#305" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">305</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Created alias directories:                                                                 <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/symlinks.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">symlinks.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/symlinks.py#223" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">223</span></a>
         <span style="color: #800080; text-decoration-color: #800080">/tmp/refgenie_demo_2hzk45es/alias/Bacteriophage-T7/fasta/</span><span style="color: #ff00ff; text-decoration-color: #ff00ff">default</span>                           <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Created alias directories:                                                                 <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/symlinks.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">symlinks.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/symlinks.py#223" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">223</span></a>
         <span style="color: #800080; text-decoration-color: #800080">/tmp/refgenie_demo_2hzk45es/alias/Bacteriophage-T7/bowtie2_index/</span><span style="color: #ff00ff; text-decoration-color: #ff00ff">2.3.5</span>                     <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Added alias: myFavGenome                                                                      <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/alias.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">alias.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/alias.py#305" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">305</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Created alias directories: <span style="color: #800080; text-decoration-color: #800080">/tmp/refgenie_demo_2hzk45es/alias/myFavGenome/fasta/</span><span style="color: #ff00ff; text-decoration-color: #ff00ff">default</span>     <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/symlinks.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">symlinks.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/symlinks.py#223" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">223</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Created alias directories:                                                                 <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/symlinks.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">symlinks.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/symlinks.py#223" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">223</span></a>
         <span style="color: #800080; text-decoration-color: #800080">/tmp/refgenie_demo_2hzk45es/alias/myFavGenome/bowtie2_index/</span><span style="color: #ff00ff; text-decoration-color: #ff00ff">2.3.5</span>                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>



The new alias should be listed in the aliases:


```python
refgenie.alias.table()
```




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="font-style: italic">                             Genome aliases                             </span>
┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃<span style="font-weight: bold"> Aliases                           </span>┃<span style="font-weight: bold"> Genome digest                    </span>┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ rCRSd                             │ jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP │
│ t7, Bacteriophage-T7, myFavGenome │ kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB │
└───────────────────────────────────┴──────────────────────────────────┘
</pre>




The command not only creates a new alias, but also creates a symbolic links to the files in the `data` directory for that genome.

Conversely, alias removal will remove the symbolic links, but not the files in the `data` directory.


```python
refgenie.alias.remove("myFavGenome")
```


<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Deleting alias-owned files: <span style="color: #800080; text-decoration-color: #800080">/tmp/refgenie_demo_2hzk45es/alias/</span><span style="color: #ff00ff; text-decoration-color: #ff00ff">myFavGenome</span>                  <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/symlinks.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">symlinks.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/symlinks.py#145" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">145</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Deleting alias-owned files: <span style="color: #800080; text-decoration-color: #800080">/tmp/refgenie_demo_2hzk45es/builds/</span><span style="color: #ff00ff; text-decoration-color: #ff00ff">myFavGenome</span>                 <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/symlinks.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">symlinks.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/symlinks.py#145" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">145</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Removed alias: myFavGenome                                                                    <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/alias.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">alias.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/alias.py#337" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">337</span></a>
</pre>



## Retrieve paths to assets

Most importantly, we can retrieve paths to refgenie-managed files.

All below commands will return the same path to the fasta file managed by Refgenie:


```python
print(refgenie.asset.seek("t7", "fasta"))
print(refgenie.asset.seek("Bacteriophage-T7", "fasta", "default"))
print(refgenie.asset.seek("t7", "fasta", "default", "fasta"))
```


<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #800080; text-decoration-color: #800080">/tmp/refgenie_demo_2hzk45es/alias/t7/fasta/default/</span><span style="color: #ff00ff; text-decoration-color: #ff00ff">t7.fa</span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #800080; text-decoration-color: #800080">/tmp/refgenie_demo_2hzk45es/alias/Bacteriophage-T7/fasta/default/</span><span style="color: #ff00ff; text-decoration-color: #ff00ff">Bacteriophage-T7.fa</span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #800080; text-decoration-color: #800080">/tmp/refgenie_demo_2hzk45es/alias/t7/fasta/default/</span><span style="color: #ff00ff; text-decoration-color: #ff00ff">t7.fa</span>
</pre>



## Remove an asset

Let's remove the `bowtie2_index` asset for the `dm6` genome.


```python
t7 = refgenie.alias.resolve(GenomeAlias("t7"))
refgenie.asset.remove("bowtie2_index", "2.3.5", genome_digest=t7)

```


<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Queued asset files for removal after commit:                                                 <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/db/events.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">events.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/db/events.py#108" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">108</span></a>
         <span style="color: #800080; text-decoration-color: #800080">/tmp/refgenie_demo_2hzk45es/data/kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB/bowtie2_index/</span><span style="color: #ff00ff; text-decoration-color: #ff00ff">2b4c1900ec5a</span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">             </span>
         <span style="color: #ff00ff; text-decoration-color: #ff00ff">2b9fce39e45a71e6bed38ca234d394b8af28c2a3b0e09548d51e</span>                                         <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">             </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Removed directory:                                                                          <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/db/cleanup.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">cleanup.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/db/cleanup.py#107" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">107</span></a>
         <span style="color: #800080; text-decoration-color: #800080">/tmp/refgenie_archive_demo_3ps6xfvv/kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB/bowtie2_index/</span><span style="color: #ff00ff; text-decoration-color: #ff00ff">2.3.5</span>    <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">              </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Removed empty directory:                                                                     <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/db/events.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">events.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/db/events.py#185" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">185</span></a>
         <span style="color: #800080; text-decoration-color: #800080">/tmp/refgenie_archive_demo_3ps6xfvv/kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB/</span><span style="color: #ff00ff; text-decoration-color: #ff00ff">bowtie2_index</span>           <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">             </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Removed directory:                                                                          <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/db/cleanup.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">cleanup.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/db/cleanup.py#107" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">107</span></a>
         <span style="color: #800080; text-decoration-color: #800080">/tmp/refgenie_demo_2hzk45es/data/kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB/bowtie2_index/</span><span style="color: #ff00ff; text-decoration-color: #ff00ff">2b4c1900ec5</span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">              </span>
         <span style="color: #ff00ff; text-decoration-color: #ff00ff">a2b9fce39e45a71e6bed38ca234d394b8af28c2a3b0e09548d51e</span>                                       <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">              </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Removed asset group and all assets <span style="color: #008000; text-decoration-color: #008000">'kN9XHLKLS_u7ei2GH87H-qpQrkz8moPB/bowtie2_index'</span>         <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/manager.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">manager.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset/manager.py#568" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">568</span></a>
</pre>



## Data channels

Refgenie supports data channels, which are used to allow third-party tool developers to expose their recipes and asset classes to Refgenie ecosystem.
In the simplest case, data channels is just a github repository with an index file that lists available asset classes and recipes, like so:

```yaml
asset_class:
  dir: asset_classes # optional, needed only if the asset classes are stored in a subdirectory
  files: # list of asset class files, relative to the index file (or directory)
    - fasta.yaml 
    - bowtie2_index.yaml
recipe:
  dir: recipes # optional, needed only if the recipes are stored in a subdirectory
  files: # list of recipe files, relative to the index file (or directory)
    - fasta.yaml
    - bowtie2_index.yaml
```

One such example is the [refgenie/recipes](https://github.com/refgenie/recipes/blob/refgenie1/index.yaml) repository, which can be added as a data channel to refgenie in the following way:


```python
from refgenie.db.tables import DataChannelType

data_channel = refgenie.sources.add_channel(
    name="refgenie-recipes",
    type=DataChannelType.https,
    index_address="https://refgenie.github.io/refgenie-registry/index.yaml",
    description="Refgenie recipes channel",
)

print(refgenie.sources.channels_table())
```


<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="font-style: italic">                                                   Data Channels                                                   </span>
┏━━━━━━━━━━━━━━━━━━┳━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━┓
┃<span style="font-weight: bold"> Name             </span>┃<span style="font-weight: bold"> Type  </span>┃<span style="font-weight: bold"> Index Address                           </span>┃<span style="font-weight: bold"> Description              </span>┃<span style="font-weight: bold"> Credentials set </span>┃
┡━━━━━━━━━━━━━━━━━━╇━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━┩
│<span style="color: #008080; text-decoration-color: #008080"> refgenie-recipes </span>│ https │ https://refgenie.github.io/refgenie-re… │ Refgenie recipes channel │ False           │
└──────────────────┴───────┴─────────────────────────────────────────┴──────────────────────────┴─────────────────┘
</pre>



Subsequently, the asset classes and recipes from the data channel can be listed and added to the refgenie instance.


```python
for asset_class in refgenie.sources.iter_asset_classes("refgenie-recipes"):
    try:
        refgenie.asset_class.add(asset_class)
    except Exception as e:
        print(e)

for recipe in refgenie.sources.iter_recipes("refgenie-recipes"):
    try:
        refgenie.recipe.add(recipe)
    except Exception as e:
        print(e)
        

```


<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span> <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/index.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 </span>   <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #008000; text-decoration-color: #008000">OK"</span>                                                                                        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/abundant_sequences.yaml</span>                <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/abundant_sequences.yaml</span>         <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'abundant_sequences'</span> asset class                                              <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL: <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/bed.yaml</span>        <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span> <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/bed.yaml</span>      <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'bed'</span> asset class                                                             <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/bismark_bt2_index.yaml</span>                 <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/bismark_bt2_index.yaml</span>          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'bismark_bt2_index'</span> asset class                                               <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL: <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/blacklist.yaml</span>  <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/blacklist.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 </span>   <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">OK"</span>                                                                                        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'blacklist'</span> asset class                                                       <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/bowtie1_index.yaml</span>                     <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/bowtie1_index.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 </span>   <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">200 OK"</span>                                                                                    <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'bowtie1_index'</span> asset class                                                   <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/bowtie2_index.yaml</span>                     <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/bowtie2_index.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 </span>   <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">200 OK"</span>                                                                                    <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace">Asset class <span style="color: #008000; text-decoration-color: #008000">'bowtie2_index'</span> version <span style="color: #008000; text-decoration-color: #008000">'0.0.1'</span> already exists.
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL: <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/bwa_index.yaml</span>  <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/bwa_index.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 </span>   <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">OK"</span>                                                                                        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'bwa_index'</span> asset class                                                       <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/cellranger_reference.yaml</span>              <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/cellranger_reference.yaml</span>       <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'cellranger_reference'</span> asset class                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL: <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/dbnsfp.yaml</span>     <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span> <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/dbnsfp.yaml</span>   <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'dbnsfp'</span> asset class                                                          <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL: <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/dbsnp.yaml</span>      <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span> <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/dbsnp.yaml</span>    <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'dbsnp'</span> asset class                                                           <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL: <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/ensembl_rb.yaml</span> <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/ensembl_rb.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 </span>  <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">OK"</span>                                                                                        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'ensembl_rb'</span> asset class                                                      <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/epilog_index.yaml</span>                      <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/epilog_index.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200</span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">OK"</span>                                                                                        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'epilog_index'</span> asset class                                                    <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL: <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/fasta.yaml</span>      <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span> <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/fasta.yaml</span>    <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace">Asset class <span style="color: #008000; text-decoration-color: #008000">'fasta'</span> version <span style="color: #008000; text-decoration-color: #008000">'0.1.0'</span> already exists.
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/fasta_index.yaml</span>                       <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/fasta_index.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 </span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">OK"</span>                                                                                        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'fasta_index'</span> asset class                                                     <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/feat_annotation.yaml</span>                   <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/feat_annotation.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 </span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">200 OK"</span>                                                                                    <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'feat_annotation'</span> asset class                                                 <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL: <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/gtf.yaml</span>        <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span> <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/gtf.yaml</span>      <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'gtf'</span> asset class                                                             <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/hisat2_index.yaml</span>                      <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/hisat2_index.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200</span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">OK"</span>                                                                                        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'hisat2_index'</span> asset class                                                    <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/kallisto_index.yaml</span>                    <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/kallisto_index.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 </span>  <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">200 OK"</span>                                                                                    <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'kallisto_index'</span> asset class                                                  <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/refgene_anno.yaml</span>                      <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/refgene_anno.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200</span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">OK"</span>                                                                                        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'refgene_anno'</span> asset class                                                    <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/salmon_index.yaml</span>                      <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/salmon_index.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200</span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">OK"</span>                                                                                        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'salmon_index'</span> asset class                                                    <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/salmon_partial_sa_index.yaml</span>           <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/salmon_partial_sa_index.yaml</span>    <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'salmon_partial_sa_index'</span> asset class                                         <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/salmon_sa_index.yaml</span>                   <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/salmon_sa_index.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 </span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">200 OK"</span>                                                                                    <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'salmon_sa_index'</span> asset class                                                 <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL: <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/small_rna.yaml</span>  <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/small_rna.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 </span>   <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">OK"</span>                                                                                        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'small_rna'</span> asset class                                                       <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL: <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/star_index.yaml</span> <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/star_index.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 </span>  <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">OK"</span>                                                                                        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'star_index'</span> asset class                                                      <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/suffixerator_index.yaml</span>                <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/suffixerator_index.yaml</span>         <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'suffixerator_index'</span> asset class                                              <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/tallymer_index.yaml</span>                    <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/tallymer_index.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 </span>  <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">200 OK"</span>                                                                                    <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'tallymer_index'</span> asset class                                                  <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL: <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/tgMap.yaml</span>      <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span> <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/asset_classes/tgMap.yaml</span>    <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'tgMap'</span> asset class                                                           <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">asset_class.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/asset_class.py#90" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">90</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span> <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/index.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 </span>   <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #008000; text-decoration-color: #008000">OK"</span>                                                                                        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/abundant_sequences/recipe.yaml</span>               <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/abundant_sequences/recipe.yaml</span>        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'abundant_sequences'</span> recipe                                                       <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL: <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/bed12/recipe.yaml</span>     <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span> <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/bed12/recipe.yaml</span>   <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'bed12'</span> recipe                                                                    <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/bismark_bt2_index/recipe.yaml</span>                <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/bismark_bt2_index/recipe.yaml</span>         <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'bismark_bt2_index'</span> recipe                                                        <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL: <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/blacklist/recipe.yaml</span> <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/blacklist/recipe.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 </span>  <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">OK"</span>                                                                                        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'blacklist'</span> recipe                                                                <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/bowtie1_index/recipe.yaml</span>                    <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/bowtie1_index/recipe.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 </span>  <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">200 OK"</span>                                                                                    <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'bowtie1_index'</span> recipe                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/bowtie2_index/recipe.yaml</span>                    <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/bowtie2_index/recipe.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 </span>  <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">200 OK"</span>                                                                                    <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace">Recipe <span style="color: #008000; text-decoration-color: #008000">'bowtie2_index'</span> version <span style="color: #008000; text-decoration-color: #008000">'0.0.1'</span> already exists.
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL: <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/bwa_index/recipe.yaml</span> <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/bwa_index/recipe.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 </span>  <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">OK"</span>                                                                                        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'bwa_index'</span> recipe                                                                <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/cellranger_reference/recipe.yaml</span>             <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/cellranger_reference/recipe.yaml</span>      <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'cellranger_reference'</span> recipe                                                     <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL: <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/dbnsfp/recipe.yaml</span>    <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span> <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/dbnsfp/recipe.yaml</span>  <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'dbnsfp'</span> recipe                                                                   <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL: <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/dbsnp/recipe.yaml</span>     <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span> <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/dbsnp/recipe.yaml</span>   <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'dbsnp'</span> recipe                                                                    <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/ensembl_gtf/recipe.yaml</span>                      <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/ensembl_gtf/recipe.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200</span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">OK"</span>                                                                                        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'ensembl_gtf'</span> recipe                                                              <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/ensembl_rb/recipe.yaml</span>                       <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/ensembl_rb/recipe.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 </span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">OK"</span>                                                                                        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'ensembl_rb'</span> recipe                                                               <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/epilog_index/recipe.yaml</span>                     <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/epilog_index/recipe.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 </span>   <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">200 OK"</span>                                                                                    <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'epilog_index'</span> recipe                                                             <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL: <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/fasta/recipe.yaml</span>     <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span> <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/fasta/recipe.yaml</span>   <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace">Recipe <span style="color: #008000; text-decoration-color: #008000">'fasta'</span> version <span style="color: #008000; text-decoration-color: #008000">'0.1.0'</span> already exists.
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/fasta_index/recipe.yaml</span>                      <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/fasta_index/recipe.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200</span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">OK"</span>                                                                                        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'fasta_index'</span> recipe                                                              <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/fasta_txome/recipe.yaml</span>                      <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/fasta_txome/recipe.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200</span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">OK"</span>                                                                                        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'fasta_txome'</span> recipe                                                              <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/feat_annotation/recipe.yaml</span>                  <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/feat_annotation/recipe.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1</span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">200 OK"</span>                                                                                    <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'feat_annotation'</span> recipe                                                          <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/gencode_gtf/recipe.yaml</span>                      <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/gencode_gtf/recipe.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200</span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">OK"</span>                                                                                        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'gencode_gtf'</span> recipe                                                              <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/hisat2_index/recipe.yaml</span>                     <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/hisat2_index/recipe.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 </span>   <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">200 OK"</span>                                                                                    <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'hisat2_index'</span> recipe                                                             <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/kallisto_index/recipe.yaml</span>                   <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/kallisto_index/recipe.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 </span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">200 OK"</span>                                                                                    <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'kallisto_index'</span> recipe                                                           <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/refgene_anno/recipe.yaml</span>                     <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/refgene_anno/recipe.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 </span>   <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">200 OK"</span>                                                                                    <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'refgene_anno'</span> recipe                                                             <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/salmon_index/recipe.yaml</span>                     <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/salmon_index/recipe.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 </span>   <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">200 OK"</span>                                                                                    <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'salmon_index'</span> recipe                                                             <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/salmon_partial_sa_index/recipe.yaml</span>          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/salmon_partial_sa_index/recipe.yaml</span>   <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'salmon_partial_sa_index'</span> recipe                                                  <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/salmon_sa_index/recipe.yaml</span>                  <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/salmon_sa_index/recipe.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1</span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">200 OK"</span>                                                                                    <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'salmon_sa_index'</span> recipe                                                          <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL: <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/small_rna/recipe.yaml</span> <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/small_rna/recipe.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 </span>  <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">OK"</span>                                                                                        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'small_rna'</span> recipe                                                                <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/star_index/recipe.yaml</span>                       <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/star_index/recipe.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 </span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">OK"</span>                                                                                        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'star_index'</span> recipe                                                               <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/suffixerator_index/recipe.yaml</span>               <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/suffixerator_index/recipe.yaml</span>        <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'suffixerator_index'</span> recipe                                                       <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL:                                                                            <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/tallymer_index/recipe.yaml</span>                   <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">        </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span>                                                                          <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/tallymer_index/recipe.yaml</span> <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 </span> <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
         <span style="color: #008000; text-decoration-color: #008000">200 OK"</span>                                                                                    <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'tallymer_index'</span> recipe                                                           <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Reading YAML from URL: <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/tgMap/recipe.yaml</span>     <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">io.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/utils/io.py#31" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">31</span></a>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> HTTP Request: <span style="color: #808000; text-decoration-color: #808000; font-weight: bold">GET</span> <span style="color: #0000ff; text-decoration-color: #0000ff; text-decoration: underline">https://refgenie.github.io/refgenie-registry/recipes/tgMap/recipe.yaml</span>   <a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">_client.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/.local/lib/python3.12/site-packages/httpx/_client.py#1025" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">1025</span></a>
         <span style="color: #008000; text-decoration-color: #008000">"HTTP/1.1 200 OK"</span>                                                                          <span style="color: #7f7f7f; text-decoration-color: #7f7f7f">               </span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace"><span style="color: #000080; text-decoration-color: #000080">INFO    </span> Registered <span style="color: #008000; text-decoration-color: #008000">'tgMap'</span> recipe                                                                    <a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">recipe.py</span></a><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">:</span><a href="file:///home/nsheff/Dropbox/workspaces/refgenie/repos/refgenie1/refgenie/managers/recipe.py#140" target="_blank"><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">140</span></a>
</pre>



Alternatively, the same can be achieved by running the following CLI command:

```bash
refgenie data-channel sync refgenie-recipes --exists-ok
```

## RefgetStore sequence collection comparison

Refgenie uses RefgetStore (from gtars) for sequence collection operations. Under the hood, refgenie uses the SeqCol digests to uniquely identify genomes. You can use RefgetStore directly to compute digests and compare sequence collections.


```python
from gtars.refget import RefgetStore

store = RefgetStore.in_memory()
store.set_quiet(True)
d1, _ = store.add_sequence_collection_from_fasta(
    REFGENIE_CODE_PATH.parent / "tests/data/rCRSd.fa"
)
d2, _ = store.add_sequence_collection_from_fasta(
    REFGENIE_CODE_PATH.parent / "tests/data/rCRSd-extra.fa"
)
store.compare(d1.digest, d2.digest)
```




    {'digests': {'a': 'jthDpfNIgzM5AGJlOkRtfnky4rXMBIUP',
      'b': 'smiTbD3jP5EwF4DNWVm0c6DGKRlFHfas'},
     'attributes': {'a_only': [],
      'b_only': [],
      'a_and_b': ['lengths',
       'name_length_pairs',
       'names',
       'sequences',
       'sorted_name_length_pairs',
       'sorted_sequences']},
     'array_elements': {'a_count': {'sorted_name_length_pairs': 1,
       'sequences': 1,
       'sorted_sequences': 1,
       'lengths': 1,
       'names': 1,
       'name_length_pairs': 1},
      'b_count': {'sequences': 2,
       'names': 2,
       'sorted_name_length_pairs': 2,
       'sorted_sequences': 2,
       'name_length_pairs': 2,
       'lengths': 2},
      'a_and_b_count': {'sorted_name_length_pairs': 1,
       'sorted_sequences': 1,
       'name_length_pairs': 1,
       'names': 1,
       'sequences': 1,
       'lengths': 1},
      'a_and_b_same_order': {'lengths': None,
       'name_length_pairs': None,
       'sorted_sequences': None,
       'names': None,
       'sequences': None,
       'sorted_name_length_pairs': None}}}


