## NOTE: this snakemake file has been generated from a template. The template requires the following variables to be defined for rendeting:
# - asset_group_names: a list of asset names to be built
# - asset_build_rule_specs: a list of asset build rule specifications
#   - asset_name: the name of the asset to be built
#   - input_asset_names: a list of asset names that are required to build the asset
#   - input_files: a list of input files that are required to build the asset
#   - params: a dictionary of parameters to be passed to the refgenie build command
#   - needs_genome_init: whether the rule depends on genome initialization
# - has_genome_init_rule: whether to include the genome_init rule

from typing import Optional
from refgenie import Refgenie
from refgenie.models import BuildCommandValues


# workflow configuration
configfile: "config.yaml"
envvars:
    "REFGENIE_INPUTS" # used in PEP, handy to be enforced here
# sample configuration
pepfile: "pep/config.yaml"

print(pep)
print(pep.sample_table)

r = Refgenie()

# Asset groups the PEP actually asks us to build. Resolving a default asset name
# means RUNNING the recipe's tool (see below), so we only pay that cost -- and
# only enforce correctness -- for recipes with real build targets.
QUEUED_ASSET_GROUPS = set()
for _sample_name in pep.sample_table["sample_name"]:
    _groups = pep.get_sample(sample_name=_sample_name).asset_group_name
    if isinstance(_groups, str):
        _groups = [_groups]
    QUEUED_ASSET_GROUPS.update(_groups)

# (recipe_name, recipe_version) -> resolved name. Resolution shells out, and
# get_all_asset_targets asks for the same names again per genome.
_ASSET_NAME_CACHE = {}

def get_default_asset_name_from_recipe(
    refgenie: Refgenie,
    recipe_name: str,
    recipe_version: Optional[str] = None,
    fallback: str = "default",
):
    """Resolve a recipe's default asset name.

    Assets are named after the version of the tool that built them, and the only
    way to learn that version is to run the tool -- recipes declare it as a shell
    command in `custom_seek_keys`. Snakemake needs the name at parse time because
    it is a component of the rule's output path.

    For a QUEUED recipe an unresolvable name is fatal. The name resolved here is
    interpolated into the `refgenie build ...:<name>` shell command, and the
    builder does not re-resolve a name it was given, so falling back to
    `default` would build and push a mis-named asset to the public bucket --
    losing the version provenance the naming scheme exists to record, and
    orphaning the correctly-named asset built by any later successful run.

    For an un-queued recipe the name is never used to build anything, so we skip
    resolution entirely rather than launch a container to name a target nobody
    asked for.
    """
    cache_key = (recipe_name, recipe_version)
    if cache_key in _ASSET_NAME_CACHE:
        return _ASSET_NAME_CACHE[cache_key]

    if recipe_name not in QUEUED_ASSET_GROUPS:
        _ASSET_NAME_CACHE[cache_key] = fallback
        return fallback

    recipe = refgenie.recipe.get(recipe_name, recipe_version)
    try:
        resolved_custom_seek_keys = (
            refgenie.resolve_custom_seek_keys(recipe.custom_seek_keys or {}) or {}
        )
        namespaces = BuildCommandValues(
            custom_seek_keys=resolved_custom_seek_keys,
            asset_group_name=recipe_name,
            genome_digest="",
            genome_folder=Path(""),
        )
        default_asset_name = refgenie.resolve_default_asset(
            default_asset=recipe.default_asset, namespaces=namespaces
        )
    except Exception as e:
        raise RuntimeError(
            f"Cannot resolve the default asset name for queued recipe "
            f"'{recipe_name}': {e}. Refusing to fall back to '{fallback}' -- "
            f"that would build and publish a mis-named asset. Check that the "
            f"recipe's custom_seek_keys commands are executable here, and that "
            f"its default_asset template references a key it actually defines."
        ) from e

    if not default_asset_name:
        raise RuntimeError(
            f"Queued recipe '{recipe_name}' resolved an empty default asset "
            f"name from default_asset={recipe.default_asset!r}. Refusing to "
            f"fall back to '{fallback}' -- that would publish a mis-named asset."
        )

    _ASSET_NAME_CACHE[cache_key] = default_asset_name
    return default_asset_name

# get the default asset (previously tag) names from the recipe.
{%-for asset_name in asset_names %}
{{asset_name}}_asset_name = get_default_asset_name_from_recipe(r, "{{asset_name}}")
{%-endfor %}

# retrieve paths to the target files for the assets.
# these files are an artificial construct in refgenie and fit well with the snakemake model to build the DAG.
{%-for asset_name in asset_names %}
{{asset_name}}_target = r.get_asset_build_target_template("{{asset_name}}", {{asset_name}}_asset_name)
{%-endfor %}
{% if has_genome_init_rule %}
# genome init target — every asset build depends on this sentinel file.
genome_init_target = r.get_genome_init_target_template()
{%- endif %}

def get_all_asset_targets(wildcards):
    """
    This is a snakmake 'input function', which gets all asset
    targets to be build based on the PEP configuration.
    """
    assets = []
    for genome_name in pep.sample_table["genome_name"]:
        assets.extend(
            [
                r.get_asset_build_target_template(
                    asset_group_name, get_default_asset_name_from_recipe(r, asset_group_name)
                )
                .as_posix()
                .format(genome_name=genome_name)
                for asset_group_name in pep.get_sample(sample_name=genome_name).asset_group_name
            ]
        )
    return assets

rule all:
    input:
        get_all_asset_targets
    shell:
        "refgenie1 list" # any command will do here
{% if has_genome_init_rule %}
rule genome_init:
    input:
        # FASTA only. Genome METADATA is deliberately NOT an input here: a file in
        # `input:` is a rebuild trigger, so listing the FHR sidecar made a
        # description edit re-init the genome and cascade a full asset rebuild.
        # Metadata is applied after the build fan-out instead (`refgenie1 genome
        # set-metadata --fhr`), which touches no build outputs.
        fasta_file_path = lambda w: pep.get_sample(sample_name=w.genome_name).fasta_file_path,
    output:
        genome_init_target = genome_init_target
    resources:
        # Serializes genome_init against every other genome_init. RefgetStore takes
        # no file lock and rewrites its index files (collections.rgci and the alias
        # TSV) WHOLESALE from the in-memory state loaded when the store was opened.
        # Two concurrent inits are therefore last-writer-wins: the loser's
        # collection is left on disk under collections/ but vanishes from the index,
        # so its alias stops resolving and every asset build for it fails with
        # "Genome not found". This is silent -- both inits report success.
        # The cap only bites if the workflow is run with a matching
        # `--resources refget_store_writer=1` (see the Rivanna profile); snakemake
        # treats an undeclared custom resource as unlimited.
        refget_store_writer = 1
    shell:
        "refgenie1 genome init --fasta {input.fasta_file_path} --name {wildcards.genome_name} --force && touch {output.genome_init_target}"
{%- endif %}
{%-for asset_build_rule_spec in asset_build_rule_specs %}
rule build_{{asset_build_rule_spec.asset_name}}:
    input:
        {%-if asset_build_rule_spec.needs_genome_init %}
        genome_init_target = genome_init_target,
        {%-endif %}
        {%-for input_asset_name in asset_build_rule_spec.input_asset_names %}
        {{input_asset_name}}_target = {{input_asset_name}}_target,
        {%-endfor %}
        {%-for input_file in asset_build_rule_spec.input_files %}
        {{input_file}}_file_path = lambda w: pep.get_sample(sample_name=w.genome_name).{{input_file}}_file_path,
        {%-endfor %}
    output:
        {{asset_build_rule_spec.asset_name}}_target = {{asset_build_rule_spec.asset_name}}_target
    params:
        {%-for param_name, param_value in asset_build_rule_spec.input_params.items() %}
        {{param_name}} = '{{param_value}}',
        {%-endfor %}
    shell:
        "refgenie1 build {wildcards.genome_name}/{{asset_build_rule_spec.asset_name}}:{%raw%}{{%endraw%}{{asset_build_rule_spec.asset_name}}_asset_name} --stage {%-for input_file in asset_build_rule_spec.input_files %} --files {{input_file}}={input.{{input_file}}_file_path}{%-endfor %}{%-for param_name, param_value in asset_build_rule_spec.input_params.items() %} --param {{param_name}}={{param_value}}{%-endfor %}"
{%-endfor %}