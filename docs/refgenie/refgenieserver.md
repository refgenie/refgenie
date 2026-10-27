# Set up your own refgenie server

## Why

We don't anticipate many people wanting to run their own servers. Typically, you'll only need to run refgenie either from the command-line or via the Python API. These clients can interact with existing refgenie servers to pull down genome assets. But what if you want to build your own server? There are a few use cases where this can make sense.

First, perhaps you want a private, local server running on your internal network. This could speed up access to private refgenie assets across an organization. Or, you may want to make some particular assets available to the community. Building on the refgenie infrastructure will simplify distribution and make it so that your users can download your resource through a familiar interface.

## How

It's pretty simple: the software that runs refgenie server is [available on GitHub](http://github.com/databio/refgenieserver). There, you will find detailed instructions on how to run it yourself. In a nutshell, you'll just run a docker command like this:

```
docker run --rm -d -p 80:80 \
	-v /path/to/genomes_archive:/genomes \
	--name refgenieservercon \
	refgenieserverim refgenieserver serve -c /genomes/genome_config.yaml
```

Mount your archived genomes folder to `/genomes` in the container, and you're essentially good to go.

### References

We have scripted the process of building and archiving the assets to serve with `refgenieserver`. The process usually includes the following steps:

1. Download raw input files for assets (FASTA files, GTF files etc.)
2. Build assets with refgenie build in a local refgenie instance
3. Archive assets with refgenieserver archive
4. Deploy the server (run `refgenieserver serve` on a cluster or locally)

Check out these GitHub repositories for more details and all the required metadata:

- [`refgenie/refgenomes.databio.org`](https://github.com/refgenie/refgenomes.databio.org) (primary refgenie assets server instance)
- [`refgenie/rg.databio.org`](https://github.com/refgenie/rg.databio.org) (development refgenie assets server instace)
