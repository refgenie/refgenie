import { useRecipes } from '../hooks/queries/useRecipes';
import { ActionBar } from '../components/common/ActionBar';
import { ScopedLink } from '../components/common/ScopedLink';
import { SearchListPage } from '../components/common/SearchListPage';
import { RECIPE_SEARCH_FIELDS } from '../services/resources/recipes';
import type { Column } from '../components/common/DataTable';
import type { RecipeResponse } from '../types/api';

const NA = <span className="rg-muted">NA</span>;

const COLUMNS: Array<Column<RecipeResponse>> = [
  {
    key: 'name',
    header: 'Name',
    render: (recipe) =>
      recipe.id !== null ? (
        <ScopedLink className="rg-link font-medium" to={`/recipes/${recipe.id}`}>
          {recipe.name}
        </ScopedLink>
      ) : (
        <span>{recipe.name}</span>
      ),
  },
  { key: 'version', header: 'Version', render: (recipe) => recipe.version },
  {
    key: 'description',
    header: 'Description',
    className: 'rg-table__cell--wrap',
    render: (recipe) => recipe.description ?? NA,
  },
  {
    key: 'default_asset',
    header: 'Default asset',
    className: 'rg-table__cell--wrap',
    render: (recipe) => recipe.default_asset,
  },
  {
    key: 'docker_image',
    header: 'Docker image',
    className: 'rg-table__cell--wrap',
    render: (recipe) => recipe.docker_image ?? NA,
  },
];

export function RecipesPage() {
  return (
    <SearchListPage
      title="Recipes"
      actions={<ActionBar slot="recipe" />}
      lede={
        <>
          A recipe is the versioned build script for one asset class: the inputs it needs, the
          commands it runs, and the container it runs them in. It is what turns a genome plus
          some input files into a finished asset.
        </>
      }
      useList={useRecipes}
      searchFields={RECIPE_SEARCH_FIELDS}
      columns={COLUMNS}
      rowKey={(recipe, index) => String(recipe.id ?? `${recipe.name}-${index}`)}
      emptyMessage="No recipes registered."
    />
  );
}
