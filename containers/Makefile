refgenie:
	mkdir -p includes
	wget https://github.com/nsheff/docker/raw/master/includes/twoBitToFa -P includes
	wget https://github.com/nsheff/docker/blob/master/includes/epilog_indexer.py -P includes
	docker build -t databio/refgenie -f Dockerfile_refgenie .
