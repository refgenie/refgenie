"""
Constant variables for refgenie package.
Ones that are integral to refgenconf and/or refgenieserver should be defined in refgenconf.const
"""
from refgenconf.const import *

BUILD_CMD = "build"
INIT_CMD = "init"
PULL_CMD = "pull"
LIST_LOCAL_CMD = "list"
LIST_REMOTE_CMD = "listr"
GET_ASSET_CMD = "seek"
INSERT_CMD = "add"
REMOVE_CMD = "remove"
GETSEQ_CMD = "getseq"
TAG_CMD = "tag"
ID_CMD = "id"
SUBSCRIBE_CMD = "subscribe"
UNSUBSCRIBE_CMD = "unsubscribe"

GENOME_ONLY_REQUIRED = [REMOVE_CMD, GETSEQ_CMD]

# For each asset we assume a genome is also required
ASSET_REQUIRED = [PULL_CMD, GET_ASSET_CMD, BUILD_CMD, INSERT_CMD, TAG_CMD, ID_CMD]

SUBPARSER_MESSAGES = {
    INIT_CMD: "Initialize a genome configuration.",
    LIST_LOCAL_CMD: "List available local assets.",
    LIST_REMOTE_CMD: "List available remote assets.",
    PULL_CMD: "Download assets.",
    BUILD_CMD: "Build genome assets.",
    GET_ASSET_CMD: "Get the path to a local asset.",
    INSERT_CMD: "Add local asset to the config file.",
    REMOVE_CMD: "Remove a local asset.",
    GETSEQ_CMD: "Get sequences from a genome.",
    TAG_CMD: "Tag an asset.",
    ID_CMD: "Return the asset digest.",
    SUBSCRIBE_CMD: "Add a refgenieserver URL to the config.",
    UNSUBSCRIBE_CMD: "Remove a refgenieserver URL from the config."
}
