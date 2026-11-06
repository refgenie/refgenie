import { listQuery } from '../../utils/pagination';
import type { ApiClient, RequestInitLite } from '../http';
import type { RecipeResponse } from '../../types/api';
import type { ListParams, PageResponse } from '../../types/pagination';

export const RECIPE_SEARCH_FIELDS = ['name', 'version', 'description'] as const;

/** One page big enough for every recipe on a server; the API caps at 1000. */
export const RECIPE_LIST_LIMIT = 200;

export interface ListRecipesParams extends ListParams {
  name?: string;
  version?: string;
  /** Filters on `AssetClass.name`, not id. */
  output_asset_class?: string;
}

export const listRecipes = (
  c: ApiClient,
  p: ListRecipesParams = {},
  init?: RequestInitLite,
) =>
  c.get<PageResponse<RecipeResponse>>(
    '/recipes',
    {
      name: p.name,
      version: p.version,
      output_asset_class: p.output_asset_class,
      ...listQuery(p),
    },
    init,
  );

export const getRecipe = (c: ApiClient, id: number, init?: RequestInitLite) =>
  c.get<RecipeResponse>(`/recipes/${id}`, undefined, init);
