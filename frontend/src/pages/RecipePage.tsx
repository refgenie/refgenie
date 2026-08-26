import { useParams } from 'react-router-dom';
import { useRecipe } from '../hooks/queries/useRecipes';
import { useUiConfig } from '../hooks/useUiConfig';
import { Breadcrumbs } from '../components/common/Breadcrumbs';
import { ErrorState, LoadingState } from '../components/common/states';
import { RecipeDetailPanel } from '../components/recipes/RecipeDetailPanel';

export function RecipePage() {
  const { id } = useParams<{ id: string }>();
  const config = useUiConfig();
  const numericId = Number(id);
  const query = useRecipe(Number.isFinite(numericId) ? numericId : undefined);

  if (query.isPending) return <LoadingState label="Loading recipe" />;
  if (query.error) {
    return <ErrorState error={query.error} subject={id} apiBase={config.api_base} />;
  }
  if (!query.data) return null;

  return (
    <div className="flex flex-col gap-8">
      <Breadcrumbs items={[{ label: 'Recipes', to: '/recipes' }, { label: query.data.name }]} />
      <h1 className="text-3xl font-bold">{query.data.name}</h1>
      <RecipeDetailPanel recipe={query.data} />
    </div>
  );
}
