# Installing refgenie

It's possible to download [pre-built refgenie assemblies](download.md) manually, but it's easier to install and use the `refgenie` command-line interface. Install refgenie from [GitHub releases](https://github.com/databio/refgenie/releases) or from PyPI with `pip`:


```console
pip install --user refgenie
```

Update with:

```console
pip install --user --upgrade refgenie
```

See if your install worked by invoking `refgenie` from the command line:

```console
refgenie -h
```

If the `refgenie` executable in not in your `$PATH`, append this to your `.bashrc` or `.profile` (or `.bash_profile` on MACOS):

```console
export PATH=~/.local/bin:$PATH
```
