import { useQuery } from '@tanstack/react-query';
import { useApiClient } from '../useApiClient';
import { qk } from '../../services/queryKeys';
import { getRecipe, listRecipes } from '../../services/resources/recipes';
import type { ListRecipesParams } from '../../services/resources/recipes';

export function useRecipes(p: ListRecipesParams) {
  const client = useApiClient();
  return useQuery({
    queryKey: qk.recipes(p),
    queryFn: ({ signal }) => listRecipes(client, p, { signal }),
  });
}

export function useRecipe(id: number | undefined) {
  const client = useApiClient();
  return useQuery({
    queryKey: qk.recipe(id ?? -1),
    queryFn: ({ signal }) => getRecipe(client, id as number, { signal }),
    enabled: id !== undefined && Number.isFinite(id),
  });
}
