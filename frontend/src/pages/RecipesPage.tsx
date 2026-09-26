import { useSearchParamsState } from '../hooks/useSearchParamsState';
import { useRecipes } from '../hooks/queries/useRecipes';
import { MiniHero } from '../components/layout/MiniHero';
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
      <MiniHero
        title="Recipes"
        actions={<ActionBar slot="recipe" />}
        lede={
          <>
            A recipe is the versioned build script for one asset class: the inputs it needs, the
            commands it runs, and the container it runs them in. It is what turns a genome plus
            some input files into a finished asset.
          </>
        }
      />

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
