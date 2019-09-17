# FAQ

## Can I use `refgenie` with my own genome resources I've already set up?

Yes, you can. Of course, one of refgenie's strengths is that it makes it easy to start a new genome resource folder from scratch. But if you've already set yours up and want to use *other* parts of the refgenie system (like the Python API, for instance), you can also do that. All you need to do is write your assets into your genome configuration file manually. Refgenie's automatic management of the configuration file won't alter these so you can use them for whatever you need.

## Can I add an asset to refgenie server?

Not to the central server -- at least, not automatically. But what you *can* do is run your own refgenie server... Since we made the `refgenieserver` software open-source and containerized, it's simple to host your own assets either on a local or publicly accessible web server. It's also possible that if the asset you're trying to host is of broad enough appeal, we'd be willing to add it to the central server, just drop us a line. In the future we'll also have an option for you to just add a recipe to a database of recipes.

## Can I maintain multiple versions of an asset?

Yes, you can. In `refgenie v0.7.0` we've introduced [tagging](tag.md), to facilitate just that!