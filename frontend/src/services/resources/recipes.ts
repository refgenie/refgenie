import type { ApiClient, RequestInitLite } from '../http';
import type { RecipePublic } from '../../types/api';
import type { ListParams, Paginated } from '../../types/pagination';

export const RECIPE_SEARCH_FIELDS = ['name', 'version', 'description'] as const;

export interface ListRecipesParams extends ListParams {
  name?: string;
  version?: string;
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
