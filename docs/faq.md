# FAQ

## Why isn't the refgenie executable available on PATH?

By default, Python command-line programs are installed to ~/.local/bin. You can add that location to your path by appending it (export PATH=$PATH:~/.local/bin). Add this to your `.bashrc` if you want it to persist.

## Can I use `refgenie` with my own genome resources I've already set up?

Yes, you can. Of course, one of refgenie's strengths is that it makes it easy to start a new genome resource folder from scratch. But if you've already set yours up and want to use *other* parts of the refgenie system (like the Python API, for instance), you can also do that. All you need to do is write your assets into your genome configuration file, which is easy using [refgenie add](custom_assets.md).

## Can I add an asset to refgenie server?

Not to the central server -- at least, not automatically. But what you *can* do is run your own refgenie server... Since we made the `refgenieserver` software open-source and containerized, it's simple to host your own assets either on a local or publicly accessible web server. It's also possible that if the asset you're trying to host is of broad enough appeal, we'd be willing to add it to the central server, just drop us a line. In the future we'll also have an option for you to just add a recipe to a database of recipes.

## Can I maintain multiple versions of an asset?

Yes, you can. In `refgenie v0.7.0` we've introduced [tagging](tag.md), to facilitate just that!

## Can multiple users share a single refgenie configuration file?

Yes, this is now the recommended way to use refgenie for groups. Starting with release `v0.7.0`, refgenie now supports genome config file locks and race-free writes, so refgenie will now automatically control multi-user conflicts to prevent metadata loss. With this change, multiple users can simultaneously read and write a single group-level configuration file.

## How can I track how a downloaded asset was created?

Starting with the server API `v2`, you can use an endpoint that will provide a detailed log output: `/v2/asset/{genome}/{asset}/log`. This log file specifies exactly how the asset was created.
