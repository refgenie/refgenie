/**
 * The two joins the recipe API does not do for us.
 *
 * `RecipePublic` carries only `output_asset_class_id`, so the class NAME is
 * resolved client-side against `/v4/asset_classes`; and the recipe list comes
 * back unordered, so multiple versions of one recipe are grouped and sorted
 * here rather than being presented as unrelated options.
 */

import { useMemo } from 'react';
import { useAssetClasses } from '../../hooks/queries/useAssetClasses';
import { useRecipes } from '../../hooks/queries/useRecipes';
import type { RecipePublic } from '../../types/api';

export const RECIPE_LIST_LIMIT = 200;

/** Descending semver, falling back to a string compare for non-semver tags. */
export function compareVersionsDesc(a: string, b: string): number {
  const parse = (value: string) => value.split('.').map((part) => Number.parseInt(part, 10));
  const left = parse(a);
  const right = parse(b);
  const usable = left.every(Number.isFinite) && right.every(Number.isFinite);
  if (!usable) return b.localeCompare(a);
  for (let i = 0; i < Math.max(left.length, right.length); i += 1) {
    const diff = (right[i] ?? 0) - (left[i] ?? 0);
    if (diff !== 0) return diff;
  }
  return 0;
}

export interface RecipeCatalog {
  /** Recipe name -> its versions, newest first. */
  byName: Map<string, RecipePublic[]>;
  isPending: boolean;
  error: unknown;
  outputClassName: (recipe: RecipePublic) => string;
}

export function useRecipeCatalog(): RecipeCatalog {
  const recipes = useRecipes({ limit: RECIPE_LIST_LIMIT });
  const assetClasses = useAssetClasses({ limit: RECIPE_LIST_LIMIT });

  const recipeItems = recipes.data?.items;
  const classItems = assetClasses.data?.items;

  return useMemo(() => {
    const classNames = new Map<number, string>();
    for (const cls of classItems ?? []) {
      if (cls.id !== null && cls.id !== undefined) classNames.set(cls.id, cls.name);
    }

    const byName = new Map<string, RecipePublic[]>();
    for (const recipe of recipeItems ?? []) {
      const list = byName.get(recipe.name) ?? [];
      list.push(recipe);
      byName.set(recipe.name, list);
    }
    for (const list of byName.values()) {
      list.sort((a, b) => compareVersionsDesc(a.version, b.version));
    }

    return {
      byName,
      isPending: recipes.isPending || assetClasses.isPending,
      error: recipes.error ?? assetClasses.error,
      outputClassName: (recipe: RecipePublic) =>
        classNames.get(recipe.output_asset_class_id) ??
        `class #${recipe.output_asset_class_id}`,
    };
  }, [
    recipeItems,
    classItems,
    recipes.isPending,
    recipes.error,
    assetClasses.isPending,
    assetClasses.error,
  ]);
}
