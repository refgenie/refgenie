import type { ApiClient, RequestInitLite } from '../http';
import type { RecipePublic } from '../../types/api';
import type { ListParams, Paginated } from '../../types/pagination';

export const RECIPE_SEARCH_FIELDS = ['name', 'version', 'description'] as const;

/** One page big enough for every recipe on a server; the API caps at 1000. */
export const RECIPE_LIST_LIMIT = 200;

export interface ListRecipesParams extends ListParams {
  name?: string;
  version?: string;
  /** Filters on `AssetClass.name`, not id — see the guard in AssetClassPage. */
  output_asset_class?: string;
}

export const listRecipes = (
  c: ApiClient,
  p: ListRecipesParams = {},
  init?: RequestInitLite,
) =>
  c.get<Paginated<RecipePublic>>(
    '/recipes',
    {
      name: p.name,
      version: p.version,
      output_asset_class: p.output_asset_class,
      q: p.q,
      search_fields: p.searchFields,
      operator: p.operator,
      offset: p.offset,
      limit: p.limit,
    },
    init,
  );

export const getRecipe = (c: ApiClient, id: number, init?: RequestInitLite) =>
  c.get<RecipePublic>(`/recipes/${id}`, undefined, init);
