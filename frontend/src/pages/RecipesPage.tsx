import { useSearchParamsState } from '../hooks/useSearchParamsState';
import { useRecipes } from '../hooks/queries/useRecipes';
import { RecipeTable } from '../components/recipes/RecipeTable';
import { SearchBox } from '../components/common/SearchBox';
import { Pagination } from '../components/common/Pagination';
import { ActionBar } from '../components/common/ActionBar';
import { EmptyState } from '../components/common/states';
import { RECIPE_SEARCH_FIELDS } from '../services/resources/recipes';

export function RecipesPage() {
  const [state, actions] = useSearchParamsState();
  const query = useRecipes({
    q: state.q || undefined,
    searchFields: state.fields.length ? state.fields : undefined,
    operator: state.operator,
    offset: state.offset,
    limit: state.limit,
  });

  return (
    <div className="flex flex-col gap-6">
      <header className="flex items-start justify-between gap-4 flex-wrap">
        <h1 className="text-3xl font-bold">Recipes</h1>
        <ActionBar slot="recipe" />
      </header>

      <SearchBox
        value={state.q}
        onChange={actions.setQuery}
        fields={RECIPE_SEARCH_FIELDS}
        selectedFields={state.fields}
        onFieldsChange={actions.setFields}
        operator={state.operator}
        onOperatorChange={actions.setOperator}
        placeholder="Search recipes…"
        label="Search recipes"
      />

      <RecipeTable
        recipes={query.data?.items}
        loading={query.isPending}
        error={query.error}
        onRetry={() => query.refetch()}
        empty={
          <EmptyState
            query={state.q || undefined}
            message="No recipes registered."
            onClearSearch={actions.clearSearch}
          />
        }
      />

      <Pagination pagination={query.data?.pagination} onOffsetChange={actions.setOffset} />
    </div>
  );
}
