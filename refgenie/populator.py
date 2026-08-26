"""Looper integration surface for refgenie.

:func:`looper_refgenie_populate_local` is a pre-submit populator hook: given a
``namespaces`` mapping, it returns a dict mapping every local genome to its
asset paths, so pipeline-interface templates can reference
``{ refgenie[sample.genome].fasta.chrom_sizes }``.

The registry-path resolution this builds on lives in
:mod:`refgenie.core.populate`; this module is the thin wrapper users name by
its dotted path (``refgenie.populator.looper_refgenie_populate_local``) in
their own configuration. Nothing here imports a workflow engine, so the
package carries no dependency on any submission tool.
"""

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from refgenie.core import Refgenie


def looper_refgenie_populate_local(namespaces: Mapping[str, Any]) -> dict[str, Any]:
    """Build ``namespaces['refgenie']`` from a local refgenie1 database.

    Pre-submit populator hook mirroring the contract of the legacy
    ``refgenconf`` populator: walks every local genome and every
    (asset_group, asset, seek_key) leaf, and writes a flat namespace
    dict shaped as ``{genome_alias: {asset_group: {seek_key: path_str}}}``
    into ``namespaces['refgenie']``. Pipeline-interface Jinja templates
    can then reference ``refgenie[sample.genome].<asset_group>.<seek_key>``
    directly.

    For each (genome, asset_group), the *default* asset is selected via
    :py:meth:`refgenie.managers.asset.manager.AssetManager.get_default`
    and that asset's seek_keys are emitted. Non-default assets are not
    walked — the legacy populator's tag selection collapses to default
    in refgenie1 because there is no per-asset tag concept.

    Connection resolution order:

    1. ``namespaces['pipeline']['var_templates']['refgenie_db_config']``
       if set and non-empty (path to ``refgenie_db_config.yaml``).
    2. ``$REFGENIE_DB_CONFIG_PATH`` environment variable (refgenie1's
       own default; honored implicitly by the ``Refgenie()`` constructor).

    Args:
        namespaces: Namespaces mapping supplied by the caller. A
            ``pipeline`` key (with optional ``var_templates``) is
            consulted if present. Mutated in place (``refgenie`` and
            ``pipeline`` keys are set/updated) in addition to being
            returned — see the merge-contract comment below for why.

    Returns:
        dict: Mapping with ``refgenie`` (the populated namespace) and
        any pre-existing ``refgenie://`` URI strings resolved via
        :py:meth:`Refgenie.populate`.

    Notes:
        - Returns ``str`` (not ``Path``) for every leaf — refgenie1's
          ``asset.seek`` returns ``Path``, which can confuse downstream
          callers expecting refgenconf semantics.
        - Failure to resolve a single seek_key (missing asset, missing
          file) is logged but does not abort the whole populator: the
          missing leaf is simply omitted, mirroring refgenconf behavior
          where Jinja ``{% if refgenie[g].x is defined %}`` guards are
          the contract for opt-in assets.
        - ``project.refgenie.tag_overrides`` and ``path_overrides`` from
          the legacy populator are NOT honored here. Refgenie1 has no
          per-asset tag concept; ``tag_overrides`` is meaningless. If
          required, file an issue and the contract can be revisited.
    """
    pipeline = namespaces.get("pipeline", {}) or {}
    var_templates = pipeline.get("var_templates", {}) or {}
    db_config_path = var_templates.get("refgenie_db_config") or None
    # Empty string from env-var substitution counts as "not set"
    if isinstance(db_config_path, str) and not db_config_path.strip():
        db_config_path = None
    # Refgenie() expects Path | None — coerce to Path so a bare str
    # from a pipeline-interface var_templates lookup works.
    if isinstance(db_config_path, str):
        db_config_path = Path(db_config_path)

    r = Refgenie(database_config_path=db_config_path)

    paths_dict: dict[str, dict[str, dict[str, str]]] = {}

    # Walk every local genome alias. r.alias.list_all() returns Alias rows
    # (with .name and .genome.digest). One alias may share a digest with
    # other aliases on the same genome; we emit the namespace under each
    # alias name so Jinja `refgenie[sample.genome]` works whichever name
    # the user supplies.
    for alias in r.alias.list_all():
        genome_alias: str = alias.name
        if not genome_alias:
            continue
        per_genome: dict[str, dict[str, str]] = {}

        # List every asset group on this genome
        try:
            groups = list(r.asset.list_groups(genome_names=[genome_alias]))
        except Exception:
            groups = []

        for group in groups:
            group_name = getattr(group, "name", None)
            if not group_name:
                continue

            # Pick the default asset for this group
            try:
                default_asset_name = r.asset.get_default(
                    asset_group_name=group_name, genome_name=genome_alias
                )
            except Exception:
                default_asset_name = "default"

            # Resolve the default asset through the name table (get() joins
            # assetname), rather than string-matching Asset.name, which is only
            # the publication name and need not equal the default name.
            asset_obj = None
            try:
                asset_obj = r.asset.get(
                    genome_name=genome_alias,
                    asset_group_name=group_name,
                    asset_name=default_asset_name,
                )
            except Exception:
                asset_obj = None
            if asset_obj is None:
                # Fall back to any asset in the group.
                try:
                    assets = list(
                        r.asset.list_assets(
                            genome_names=[genome_alias],
                            asset_group_name=group_name,
                        )
                    )
                except Exception:
                    assets = []
                if assets:
                    asset_obj = assets[0]
                    default_asset_name = getattr(asset_obj, "name", "default")
            if asset_obj is None:
                continue

            seek_keys_iter = getattr(asset_obj, "seek_keys", []) or []
            per_seek: dict[str, str] = {}
            for sk in seek_keys_iter:
                sk_name = getattr(sk, "name", None)
                if not sk_name:
                    continue
                try:
                    val = r.asset.seek(
                        genome_alias,
                        group_name,
                        asset_name=default_asset_name,
                        seek_key_name=sk_name,
                    )
                except Exception:
                    # Missing file / unresolvable seek_key — skip the leaf
                    continue
                per_seek[sk_name] = str(val)
            if per_seek:
                per_genome[group_name] = per_seek

        if per_genome:
            paths_dict[genome_alias] = per_genome

    # Callers merge the returned dict one leaf at a time, doing
    # x[namespace][key] = val — so x[namespace] must already exist on
    # the input namespaces before that merge. We mutate the input dict
    # to set the 'refgenie' namespace (matching the legacy populator's
    # behavior), then also include it in the returned dict so the
    # per-key merge path works correctly.
    if isinstance(namespaces, dict):
        namespaces["refgenie"] = paths_dict
    out: dict[str, Any] = {"refgenie": paths_dict}

    # Refgenconf populator also resolves embedded ``refgenie://`` URIs in
    # the namespaces. Refgenie1 has Refgenie.populate for that exact
    # purpose. We resolve only the pipeline namespace because that's
    # where pipeline-interface authors might embed registry paths.
    pipeline_block = namespaces.get("pipeline")
    if pipeline_block is not None:
        try:
            populated_pipeline = r.populate(dict(pipeline_block))
            if isinstance(populated_pipeline, dict) and isinstance(
                namespaces.get("pipeline"), dict
            ):
                # Merge populated values back into the in-place pipeline
                # namespace so the caller's update path can see them.
                namespaces["pipeline"].update(populated_pipeline)
                out["pipeline"] = populated_pipeline
        except Exception:
            # If population fails, leave the pipeline namespace alone
            pass

    return out
