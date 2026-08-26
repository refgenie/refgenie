/**
 * Recipe picker.
 *
 * Two joins the backend does not do for us: `RecipePublic` carries only
 * `output_asset_class_id`, so the class name is resolved client-side; and the
 * recipe list is unordered, so multiple versions of one recipe are collapsed
 * into a single option with a version dropdown defaulting to the highest
 * semver, sorted here.
 */

import { useMemo } from 'react';
import { FormField } from '../common/FormField';
import { Badge } from '../common/Badge';
import { useRecipeCatalog } from './recipeCatalog';
import type { RecipePublic } from '../../types/api';

export interface RecipeSelectProps {
  recipeName: string | null;
  recipeVersion: string | null;
  onChange: (recipe: RecipePublic | null) => void;
  error?: string;
}

export function RecipeSelect({
  recipeName,
  recipeVersion,
  onChange,
  error,
}: RecipeSelectProps) {
  const catalog = useRecipeCatalog();

  const grouped = useMemo(() => {
    const groups = new Map<string, RecipePublic[]>();
    for (const versions of catalog.byName.values()) {
      const newest = versions[0];
      if (!newest) continue;
      const className = catalog.outputClassName(newest);
      const list = groups.get(className) ?? [];
      list.push(newest);
      groups.set(className, list);
    }
    return [...groups.entries()].sort((a, b) => a[0].localeCompare(b[0]));
  }, [catalog]);

  const versions = recipeName ? (catalog.byName.get(recipeName) ?? []) : [];
  const selected =
    versions.find((recipe) => recipe.version === recipeVersion) ?? versions[0] ?? null;

  return (
    <div className="flex flex-col gap-2">
      <FormField
        htmlFor="build-recipe"
        label="Recipe"
        required
        error={error}
        hint={
          selected?.description ?? 'Recipes are grouped by the asset class they produce.'
        }
      >
        <select
          id="build-recipe"
          className="rg-field__input"
          value={recipeName ?? ''}
          onChange={(event) => {
            const name = event.target.value;
            const list = catalog.byName.get(name);
            onChange(list?.[0] ?? null);
          }}
        >
          <option value="">
            {catalog.isPending ? 'Loading recipes…' : 'Choose a recipe…'}
          </option>
          {grouped.map(([className, recipes]) => (
            <optgroup key={className} label={className}>
              {recipes.map((recipe) => (
                <option key={recipe.name} value={recipe.name}>
                  {recipe.name}
                </option>
              ))}
            </optgroup>
          ))}
        </select>
      </FormField>

      {versions.length > 1 && (
        <FormField
          htmlFor="build-recipe-version"
          label="Recipe version"
          hint="Defaults to the newest installed version."
        >
          <select
            id="build-recipe-version"
            className="rg-field__input"
            value={selected?.version ?? ''}
            onChange={(event) => {
              const match = versions.find((recipe) => recipe.version === event.target.value);
              onChange(match ?? null);
            }}
          >
            {versions.map((recipe) => (
              <option key={recipe.version} value={recipe.version}>
                {recipe.version}
              </option>
            ))}
          </select>
        </FormField>
      )}

      {selected?.docker_image && (
        <p className="text-sm">
          <Badge variant="archive">docker</Badge>{' '}
          <code className="rg-code rg-code--inline">{selected.docker_image}</code>
        </p>
      )}
    </div>
  );
}
