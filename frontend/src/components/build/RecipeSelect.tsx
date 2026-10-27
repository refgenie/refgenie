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

  /*
   * One row per recipe, flat: every row in the list is a recipe the user can
   * choose. No headers -- a header is a row that says nothing selectable, and
   * a list of 29 recipes that spends rows on 27 classes reads as doubled.
   *
   * The asset class still earns its place, inline, but only where it differs
   * from the recipe's own name (`gtf · ensembl_gtf`, `fasta · fasta_txome`),
   * because `bowtie2_index · bowtie2_index` is noise. Sorting by class then
   * name keeps recipes that produce the same class adjacent.
   */
  const entries = useMemo(() => {
    const rows: { recipe: RecipePublic; className: string; label: string }[] = [];
    for (const versions of catalog.byName.values()) {
      const newest = versions[0];
      if (!newest) continue;
      const className = catalog.outputClassName(newest);
      rows.push({
        recipe: newest,
        className,
        label: className === newest.name ? newest.name : `${className} · ${newest.name}`,
      });
    }
    return rows.sort(
      (a, b) =>
        a.className.localeCompare(b.className) || a.recipe.name.localeCompare(b.recipe.name),
    );
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
          selected?.description ??
          'Named “asset class · recipe” where the two differ; otherwise just the recipe.'
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
          {entries.map((entry) => (
            <option key={entry.recipe.name} value={entry.recipe.name}>
              {entry.label}
            </option>
          ))}
        </select>

        {/* Inside the field, not after it: the control column is where this
            annotation lines up, at any container width. */}
        {selected?.docker_image && (
          <p className="text-sm mt-1">
            <Badge variant="archive">docker</Badge>{' '}
            <code className="rg-code rg-code--inline">{selected.docker_image}</code>
          </p>
        )}
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
    </div>
  );
}
