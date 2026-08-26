/**
 * Client-side build-form validation.
 *
 * A pure function, deliberately: the schema is dynamic and per-recipe, which is
 * the case form libraries handle worst, and a pure validator is unit-testable
 * without rendering anything.
 *
 * It never duplicates the server's judgement — path existence, recipe
 * resolution and input-asset availability are `POST /actions/build/preflight`.
 * This only catches what the browser can know: something required is blank, or
 * a name contains a character that cannot survive a registry path.
 */

import type { RecipePublic } from '../../types/api';

/** Mirrors `ASSET_REGISTRY_PATH_COMPONENT_CHAR_BLACKLIST` (`refgenie/const.py`). */
export const REGISTRY_PATH_BLACKLIST = [' ', ':', '/'] as const;

export interface BuildFormValues {
  /** The alias, which `facade.build_asset` resolves itself. */
  genome: string;
  /** Carried alongside for the job target key only. */
  genomeDigest: string | null;
  recipeName: string | null;
  recipeVersion: string | null;
  assetGroupName: string;
  assetName: string;
  assetDescription: string;
  params: Record<string, string>;
  files: Record<string, string>;
  assets: Record<string, string>;
  stage: boolean;
  pullParents: boolean;
}

/**
 * Dotted keys, matching the `field` on the server's error envelope AND the
 * request body's own paths — which is why the recipe inputs are prefixed
 * `params.`: they are nested under `params` on the wire.
 */
export type FieldErrors = Record<string, string>;

/** The dotted wire path for one recipe input, used as its error key. */
export function inputFieldKey(section: 'params' | 'files' | 'assets', name: string): string {
  return `params.${section}.${name}`;
}

export function emptyBuildForm(): BuildFormValues {
  return {
    genome: '',
    genomeDigest: null,
    recipeName: null,
    recipeVersion: null,
    assetGroupName: '',
    assetName: '',
    assetDescription: '',
    params: {},
    files: {},
    assets: {},
    stage: false,
    pullParents: false,
  };
}

function blacklisted(value: string): string | null {
  const hit = REGISTRY_PATH_BLACKLIST.find((char) => value.includes(char));
  if (!hit) return null;
  return hit === ' '
    ? 'Spaces are not allowed in a registry path component.'
    : `The character "${hit}" is not allowed in a registry path component.`;
}

function entries(source: Record<string, Record<string, unknown>> | null | undefined) {
  return Object.entries(source ?? {});
}

function hasDefault(spec: Record<string, unknown>): boolean {
  const value = spec.default;
  return value !== undefined && value !== null && value !== '';
}

export function validateBuildForm(
  recipe: RecipePublic | null | undefined,
  values: BuildFormValues,
): FieldErrors {
  const errors: FieldErrors = {};

  if (!values.genome.trim()) errors.genome = 'Choose a genome.';
  if (!recipe) errors.recipe = 'Choose a recipe.';

  const group = values.assetGroupName.trim();
  if (!group) {
    errors.asset_group = 'An asset group name is required.';
  } else {
    const bad = blacklisted(group);
    if (bad) errors.asset_group = bad;
  }

  if (values.assetName) {
    const bad = blacklisted(values.assetName);
    if (bad) errors.asset = bad;
  }

  if (!recipe) return errors;

  for (const [name, spec] of entries(recipe.input_params)) {
    if (!hasDefault(spec) && !values.params[name]?.trim()) {
      errors[inputFieldKey('params', name)] = 'This parameter is required.';
    }
  }
  for (const [name, spec] of entries(recipe.input_files)) {
    if (!hasDefault(spec) && !values.files[name]?.trim()) {
      errors[inputFieldKey('files', name)] = 'A file path is required.';
    }
  }
  for (const [name, spec] of entries(recipe.input_assets)) {
    if (!hasDefault(spec) && !values.assets[name]?.trim()) {
      errors[inputFieldKey('assets', name)] = 'An input asset is required.';
    }
  }

  return errors;
}

/** `<genome>/<asset_group>:<asset_name or "auto">`, live under the form. */
export function registryPathPreview(values: BuildFormValues): string {
  const genome = values.genome || '<genome>';
  const group = values.assetGroupName || '<asset_group>';
  return `${genome}/${group}:${values.assetName || 'auto'}`;
}
