# Setting up your own refgenie server

*This is an advanced tutorial. Most users should not need to do this*.

Typically, you'll only need to run refgenie either from the command-line or via the Python API. These clients can interact with the refgenie server to pull down genome assets. But what if you want to build your own server?

Though we don't anticipate many people wanting to run their own servers, there are a few use cases where this can make sense. First, perhaps you want a private, local server running on your internal network. This could speed up access to refgenie assets. Another reason is that you may want to make some particular assets available to the community. Building on the refgenie infrastructure will simplify distribution for your, and make it so that your users can download your resource through a familiar interface. 

The software that runs refgenie server is [available on GitHub](http://github.com/databio/refgenieserver). There, you will find detailed instructions on how to run it yourself.
