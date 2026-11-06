import { ScopedLink } from '../common/ScopedLink';
import { useRecipes } from '../../hooks/queries/useRecipes';
import { RECIPE_LIST_LIMIT } from '../../services/resources/recipes';
import { ErrorState } from '../common/states';

export interface AssetClassRecipeListProps {
  /** The route's numeric id, used to drop other versions of the same name. */
  assetClassId: number;
  /** The only thing `/recipes` can filter on. */
  assetClassName: string;
}

/** The recipes that produce this asset class. */
export function AssetClassRecipeList({
  assetClassId,
  assetClassName,
}: AssetClassRecipeListProps) {
  const query = useRecipes({
    output_asset_class: assetClassName,
    pageSize: RECIPE_LIST_LIMIT,
  });

  if (query.error) return <ErrorState error={query.error} onRetry={() => query.refetch()} />;

  // The endpoint filters by class name only, so drop rows from other versions.
  const recipes = (query.data?.results ?? []).filter(
    (recipe) => recipe.output_asset_class_id === assetClassId,
  );

  return (
    <section>
      <h2 className="text-xl font-semibold mb-4">Built by</h2>
      {recipes.length === 0 ? (
        <p className="rg-muted text-sm">No recipe on this server produces this asset class.</p>
      ) : (
        <ul className="flex flex-col gap-1">
          {recipes.map((recipe, index) => (
            <li className="text-sm" key={recipe.id ?? `${recipe.name}-${index}`}>
              {recipe.id !== null ? (
                <ScopedLink className="rg-link font-medium" to={`/recipes/${recipe.id}`}>
                  {recipe.name}
                </ScopedLink>
              ) : (
                <span className="font-medium">{recipe.name}</span>
              )}
              <span className="rg-muted"> v{recipe.version}</span>
              {recipe.description && <span className="rg-muted"> — {recipe.description}</span>}
            </li>
          ))}
        </ul>
      )}
    </section>
  );
}
