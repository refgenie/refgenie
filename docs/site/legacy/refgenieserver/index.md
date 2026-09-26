# refgenieserver

[![Build Status](https://travis-ci.org/databio/refgenieserver.svg?branch=master)](https://travis-ci.org/databio/refgenieserver)

# refgenieserver

This folder contains code for an API to provide reference genomes. `refgenieserver` can do 2 things: `archive` an existing refgenie folder, and then `serve` it. 

## How to `serve`

### Building container

In the same directory as the `Dockerfile`:

```
docker build -t refgenieserverim .
```

### Running container for development:

You can run it directly after installing with `pip install`, like this:

```
refgenieserver serve -c refgenie.yaml -p 5000
```

Better, though, is to use the container. Mount a directory of files to serve at `/genomes`:

```
docker run --rm -p 80:80 --name refgenieservercon \
  -v $(pwd)/files:/genomes \
  refgenieserverim refgenieserver serve -c refgenie.yaml 
```

### Running container for production:
Run the container from the image you just built:

```
docker run --rm -d -p 80:80 \
  -v /path/to/genomes_archive:/genomes \
  --name refgenieservercon \
  refgenieserverim refgenieserver serve -c /genomes/genome_config.yaml 
```

Make sure the `genome_config.yaml` filename matches what you've named your configuration file! We use `-d` to detach so it's in background. You shouldn't need to mount the app (`-v /path/to/refgenieserver:/app`) because in this case we're running it directly. Terminate container when finished:

```
docker stop refgenieservercon
```


### Interacting with the API web server

Navigate to [http://localhost/](http://localhost/) to see the server in action.

You can see the automatic docs and interactive swagger openAPI interface at [http://localhost/docs](http://localhost/docs). That will also tell you all the endpoints, etc.


### Monitoring for errors

Attach to container to see debug output:

```
docker attach refgenieservercon
```

Grab errors:

```
docker events | grep -oP "(?<=die )[^ ]+"
```

View those error codes:

```
docker logs <error_code>
```

Enter an interactive shell to explore the container contents:

```
docker exec -it refgenieservercon sh
```

## How to `archive`

Refgenieserver can also archive your assets, creating the directory for asset archives needed to `serve`.

First, make sure the config has a `genome_archive_folder` key that points to the directory where you want to store the servable archives (`genome_archive_folder` is __not__ added automatically by [`refgenie init`](http://refgenie.databio.org/en/latest/usage/#refgenie-init-help)). Your first time you will need to manually add this to tell refgenieserver where to store the archives.

Secondly, if you wish to store the refgenieserver configuration file separately from the `genome_archive_folder`, specify a `genome_archive_config` key. The path that this key points to will be considered relative to the refgenie configuration file, unless it's absolute.

Then run: 
```
refgenieserver archive -c CONFIG
````
It just requires a `-c` argument or `$REFGENIE` environment variable.

This command will:
- create the `genome_archive` directory and structure that can be used to serve the assets
- create a server config file in that directory, which includes a couple of extra asset attributes, like `archive_digest` and `archive_size`. 

In case you already have some of the assets archived and just want to add a new one, use:

```
refgenieserver archive -c CONFIG GENOME/ASSET:TAG
```

In case you want to remove an unwanted archive, add an `-r` flag:

```
refgenieserver archive -c CONFIG -r GENOME/ASSET:TAG
```

## How to test the `refgenie` suite of software

The `refgenie` universe includes [`refgenie`](http://refgenie.databio.org/en/latest/), [`refgenconf`](http://refgenie.databio.org/en/latest/overview/#3-the-refgenconf-python-package), and [`refgenieserver`](https://github.com/databio/refgenieserver/).

The [`test_refgenie.sh`](https://github.com/databio/refgenieserver/blob/staging/test_refgenie.sh) script will test the integration of all three tools to ensure everything is functioning, particularly following any changes or updates to one, two, or all three tools.

Use it simply as follows:
```
/path/to/test_refgenie.sh
```

The script also requires Python's [virtual environment module](https://docs.python.org/3/tutorial/venv.html), [Docker](https://www.docker.com/), and [Bulker](https://bulker.databio.org/en/latest/) to successfully test all components.
