import { useParams } from 'react-router-dom';
import { useRecipe } from '../hooks/queries/useRecipes';
import { useUiConfig } from '../hooks/useUiConfig';
import { Breadcrumbs } from '../components/common/Breadcrumbs';
import { MiniHero } from '../components/layout/MiniHero';
import { ErrorState, LoadingState } from '../components/common/states';
import { RecipeDetailPanel } from '../components/recipes/RecipeDetailPanel';

export function RecipePage() {
  const { id } = useParams<{ id: string }>();
  const config = useUiConfig();
  const numericId = Number(id);
  const query = useRecipe(Number.isFinite(numericId) ? numericId : undefined);

  // No early return: the head renders in every state.
  const recipe = query.data;

  return (
    <div className="flex flex-col gap-8">
      <MiniHero
        title={recipe ? recipe.name : 'Recipe'}
        documentTitle={recipe ? undefined : false}
        breadcrumbs={
          <Breadcrumbs
            items={[
              { label: 'Recipes', to: '/recipes' },
              ...(recipe ? [{ label: recipe.name }] : []),
            ]}
          />
        }
      />

      {query.isPending && <LoadingState label="Loading recipe" />}
      {query.error && <ErrorState error={query.error} subject={id} apiBase={config.api_base} />}

      {recipe && <RecipeDetailPanel recipe={recipe} />}
    </div>
  );
}
