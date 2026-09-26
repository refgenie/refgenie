import { Link } from 'react-router-dom';
import { useRecipes } from '../../hooks/queries/useRecipes';
import { RECIPE_LIST_LIMIT } from '../../services/resources/recipes';
import { ErrorState } from '../common/states';

export interface AssetClassRecipeListProps {
  /** The route's numeric id — the exact class, used to guard the name filter. */
  assetClassId: number;
  /** The only thing `/recipes` can filter on. */
  assetClassName: string;
}

/**
 * The recipes that produce this asset class.
 *
 * This is the direction the schema supports: `Recipe.output_asset_class_id` is a
 * single non-nullable FK, so a recipe has one output class and a class has many
 * recipes. The inverse ("classes that fulfill a recipe") is not a set.
 */
export function AssetClassRecipeList({
  assetClassId,
  assetClassName,
}: AssetClassRecipeListProps) {
  const query = useRecipes({
    output_asset_class: assetClassName,
    limit: RECIPE_LIST_LIMIT,
  });

  if (query.error) return <ErrorState error={query.error} onRetry={() => query.refetch()} />;

  // `/recipes?output_asset_class=` filters on AssetClass.NAME, and AssetClass is
  // unique on (name, version) — two versions can share a name. Drop rows that
  // belong to a sibling version of this class.
  const recipes = (query.data?.items ?? []).filter(
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
                <Link className="rg-link font-medium" to={`/recipes/${recipe.id}`}>
                  {recipe.name}
                </Link>
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
