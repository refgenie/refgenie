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

  // The head renders in every state, so nothing may early-return past it: a
  // loading or missing recipe still gets a heading, a breadcrumb out, and (via
  // MiniHero) a correct document title. No explainer here -- `/recipes` above
  // defines the term, and the breadcrumb links straight back to it.
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
