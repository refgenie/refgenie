
## Docker

I have produced a docker image on DockerHub (nsheff/refgenie) that has all of these packages pre-installed, so you can run the complete indexer without worrying about paths and packages. Just clone this repo and run it with the `-d` flag. For example:

```
~/code/refgenie/src/refgenie.py -d --input rn6.fa --outfolder $HOME
```

