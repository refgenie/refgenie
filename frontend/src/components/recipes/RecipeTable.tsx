import { Link } from 'react-router-dom';
import { DataTable } from '../common/DataTable';
import type { Column } from '../common/DataTable';
import type { RecipePublic } from '../../types/api';

export interface RecipeTableProps {
  recipes: RecipePublic[] | undefined;
  loading?: boolean;
  error?: unknown;
  empty?: React.ReactNode;
  onRetry?: () => void;
}

export function RecipeTable({ recipes, loading, error, empty, onRetry }: RecipeTableProps) {
  const columns: Array<Column<RecipePublic>> = [
    {
      key: 'name',
      header: 'Name',
      render: (recipe) =>
        recipe.id !== null ? (
          <Link className="rg-link font-medium" to={`/recipes/${recipe.id}`}>
            {recipe.name}
          </Link>
        ) : (
          <span>{recipe.name}</span>
        ),
    },
    { key: 'version', header: 'Version', render: (recipe) => recipe.version },
    {
      key: 'description',
      header: 'Description',
      render: (recipe) => recipe.description ?? <span className="rg-muted">NA</span>,
    },
    {
      key: 'default_asset',
      header: 'Default asset',
      render: (recipe) => recipe.default_asset,
    },
    {
      key: 'docker_image',
      header: 'Docker image',
      render: (recipe) => recipe.docker_image ?? <span className="rg-muted">NA</span>,
    },
  ];

  return (
    <DataTable
      caption="Recipes"
      columns={columns}
      rows={recipes}
      rowKey={(recipe, index) => String(recipe.id ?? `${recipe.name}-${index}`)}
      loading={loading}
      error={error}
      empty={empty}
      onRetry={onRetry}
    />
  );
}
