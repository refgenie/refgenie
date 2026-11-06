/**
 * Client-side joins for the recipe list: resolve each output class name, and
 * group a recipe's versions, newest first.
 */

import { useMemo } from 'react';
import { useAssetClassIndex } from '../../hooks/queries/useAssetClasses';
import { useRecipes } from '../../hooks/queries/useRecipes';
import { RECIPE_LIST_LIMIT } from '../../services/resources/recipes';
import type { RecipeResponse } from '../../types/api';

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
  byName: Map<string, RecipeResponse[]>;
  isPending: boolean;
  error: unknown;
  outputClassName: (recipe: RecipeResponse) => string;
}

export function useRecipeCatalog(): RecipeCatalog {
  const recipes = useRecipes({ pageSize: RECIPE_LIST_LIMIT });
  const classIndex = useAssetClassIndex();

  const recipeItems = recipes.data?.results;
  const classData = classIndex.data;

  return useMemo(() => {
    const byName = new Map<string, RecipeResponse[]>();
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
      isPending: recipes.isPending || classIndex.isPending,
      error: recipes.error ?? classIndex.error,
      outputClassName: (recipe: RecipeResponse) =>
        classData?.byId.get(recipe.output_asset_class_id)?.name ??
        `class #${recipe.output_asset_class_id}`,
    };
  }, [
    recipeItems,
    classData,
    recipes.isPending,
    recipes.error,
    classIndex.isPending,
    classIndex.error,
  ]);
}
