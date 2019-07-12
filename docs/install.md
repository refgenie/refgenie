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

Refgenie will automatically use the config file in this environmental variable if it exists. Add this to your `.bashrc` or `.profile` if you want it to persist for future command-line sessions. You can always specify `-c` if you want to override the value in the `$REFGENIE` variable on an ad-hoc basis.

# Listing assets

Now you can use the `list` command to show local assets (which will be empty at first) or the `listr` command to show available remote assets:

```console
refgenie list
refgenie listr
```

# Populate some assets

Next you need to populate your genome folder with a few assets. You can either `pull` existing assets or `build` your own. Refgenie will manage them the same way. As an example, let's pull a bowtie2 index for a small genome, the human mitochondrial genome (it's called `rCRSd`, the "Revised Cambridge Reference Sequence" on our server).

```console
refgenie pull -g rCRSd -a bowtie2_index
```

You can also read more about [building refgenie assets](build.md).

# Seeking assets

Use the `seek` command to get paths to local assets you have already built or pulled. For example, the one we just pulled:

```console
refgenie pull -g rCRSd -a bowtie2_index
```

Or, more generally:

```console
refgenie seek -g GENOME -a ASSET
```

That's it! Explore the HOW-TO guides in the navigation bar for further details about what you can do with these functions.
