"""Dispatch table mapping parsed command models to their handlers.

Imported lazily by ``refgenie.cli.main`` after argument parsing, so the cost
of importing every command family is only paid when a command actually runs.
"""

from collections.abc import Callable

from refgenie.cli.commands import (
    alias,
    asset,
    asset_class,
    build,
    config,
    curate,
    data_channel,
    database,
    generate,
    genome,
    getseq,
    listing,
    lookup,
    populate,
    pull,
    push,
    recipe,
    remote,
    serve,
    servers,
    stage,
    store,
)

_DISPATCH: dict[type, Callable] = {
    database.InitModel: database.handle_init,
    database.PurgeModel: database.handle_purge,
    listing.ListModel: listing.handle_list,
    listing.ListrModel: listing.handle_listr,
    servers.SubscribeModel: servers.handle_subscribe,
    servers.UnsubscribeModel: servers.handle_unsubscribe,
    servers.CatalogExportModel: servers.handle_catalog_export,
    lookup.SeekModel: lookup.handle_seek,
    lookup.SeekrModel: lookup.handle_seekr,
    lookup.IdModel: lookup.handle_id,
    lookup.CompareModel: lookup.handle_compare,
    curate.RemoveModel: curate.handle_remove,
    curate.RenameModel: curate.handle_rename,
    curate.InsertModel: curate.handle_add,
    getseq.GetseqModel: getseq.handle_getseq,
    pull.PullModel: pull.handle_pull,
    pull.MirrorModel: pull.handle_mirror,
    populate.PopulateModel: populate.handle_populate,
    populate.PopulaterModel: populate.handle_populater,
    push.PushModel: push.handle_push,
    build.BuildModel: build.handle_build,
    serve.ServeModel: serve.handle_serve,
    serve.DashModel: serve.handle_dash,
    # Nested groups
    alias.AliasParser: alias.handle_alias_group,
    config.ConfigParser: config.handle_config_group,
    recipe.RecipeParser: recipe.handle_recipe_group,
    asset_class.AssetClassParser: asset_class.handle_asset_class_group,
    stage.StageParser: stage.handle_stage_group,
    data_channel.DataChannelParser: data_channel.handle_data_channel_group,
    generate.GenerateParser: generate.handle_generate_group,
    remote.RemoteParser: remote.handle_remote_group,
    store.StoreParser: store.handle_store_group,
    genome.GenomeParser: genome.handle_genome_group,
    asset.AssetParser: asset.handle_asset_group,
}


def get_dispatch() -> dict[type, Callable]:
    """Return the model-type -> handler dispatch table."""
    return _DISPATCH
