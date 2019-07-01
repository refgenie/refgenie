# Installing refgenie

Install refgenie from [GitHub releases](https://github.com/databio/refgenie/releases) or from PyPI with `pip`:

- `pip install --user refgenie`: install into user space.
- `pip install --user --upgrade refgenie`: update in user space.
- `pip install refgenie`: install into an active virtual environment.
- `pip install --upgrade refgenie`: update in virtual environment.

See if your install worked by calling `refgenie -h` on the command line. If the `refgenie` executable in not in your `$PATH`, append this to your `.bashrc` or `.profile` (or `.bash_profile` on macOS):
```console
export PATH=~/.local/bin:$PATH
```

# Initial configuration

If you're using refgenie for the first time you'll need to initialize your ***genome folder*** and configuration file. Just select a folder where you want your genome assets to live, and then try:
```console
refgenie init -c genome_folder/genome_config.yaml
```

The `refgenie` commands all require knowing where this genome config file is. You can pass it on the command line all the time (using the `-c` parameter), but this gets old. An alternative is to set up the `$REFGENIE` environment variable like so:
```console
export REFGENIE=/path/to/genome_config.yaml
```

Add this to your `.bashrc` or `.profile` if you want it to persist for future command-line sessions. You can always specific `-c` if you want to override the value in the `$REFGENIE` variable on an ad-hoc basis.

# Listing assets

Now you can use the `list` command to show local assets (which will be empty at first) or the `listr` command to show available remote assets:
```
refgenie list
refgenie listr
```

# Seeking assets

Use the `seek` command to get paths to local assets you have already built or pulled:

```console
refgenie seek -g GENOME -a ASSET
```
