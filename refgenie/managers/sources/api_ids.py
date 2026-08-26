"""OperationIds of a remote refgenie server's OpenAPI endpoints.

These name the endpoints the client resolves through the server's published
OpenAPI document, so changes of endpoint function names OR endpoints themselves
do not influence the connection. Every consumer is client-side.
"""

CUSTOM_PFX = "custom_Id"
API_ID_ALIAS_DIGEST = "get_alias_v4_aliases__name__get"
API_ID_ARCHIVE = "download_archive_v4_archives__asset_digest__download_get"
API_ID_STAGED_ASSETS = "list_staged_assets_v4_staged_assets_get"
API_ID_ASSET_ATTRS = CUSTOM_PFX + "_asset_attrs"
API_ID_GENOME_ATTRS = "get_genome_v4_genomes__digest__get"
API_ID_DIGEST = CUSTOM_PFX + "_asset_digest"
API_ID_ASSET_FILES = "list_asset_files_v4_assets__asset_digest__files_get"
API_ID_ASSET_FILE_DOWNLOAD = "download_asset_file_v4_assets__asset_digest__files__file_path__get"
API_ID_ASSET_RELATIONSHIPS = "get_asset_relationships_v4_relationships__asset_digest__get"
