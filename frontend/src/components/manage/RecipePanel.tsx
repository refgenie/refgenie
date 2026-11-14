/**
 * Installed recipes and asset classes, read-only. Adding them over HTTP would
 * allow arbitrary file reads and SSRF, so registration stays in the CLI.
 */

import { Link } from 'react-router-dom';
import { useAssetClasses } from '../../hooks/queries/useAssetClasses';
import { useRecipes } from '../../hooks/queries/useRecipes';
import { DataTable } from '../common/DataTable';
import { Badge } from '../common/Badge';
import { useRecipeCatalog } from '../build/recipeCatalog';
import type { Column } from '../common/DataTable';
import type { AssetClassPublic, RecipeResponse } from '../../types/api';

const LIMIT = 200;

function inputSummary(recipe: RecipeResponse): string {
  const counts = [
    [Object.keys(recipe.input_params ?? {}).length, 'params'],
    [Object.keys(recipe.input_files ?? {}).length, 'files'],
    [Object.keys(recipe.input_assets ?? {}).length, 'assets'],
  ] as const;
  const parts = counts.filter(([n]) => n > 0).map(([n, label]) => `${n} ${label}`);
  return parts.length ? parts.join(', ') : 'none';
}

export function RecipePanel() {
  const recipes = useRecipes({ pageSize: LIMIT });
  const catalog = useRecipeCatalog();

  const columns: Array<Column<RecipeResponse>> = [
    {
      key: 'name',
      header: 'Recipe',
      render: (recipe) =>
        recipe.id !== null && recipe.id !== undefined ? (
          <Link className="rg-link font-medium" to={`/recipes/${recipe.id}`}>
            {recipe.name}
          </Link>
        ) : (
          recipe.name
        ),
    },
    { key: 'version', header: 'Version', render: (recipe) => recipe.version },
    {
      key: 'output',
      header: 'Output class',
      render: (recipe) => catalog.outputClassName(recipe),
    },
    {
      key: 'docker',
      header: 'Docker',
      render: (recipe) =>
        recipe.docker_image ? (
          <Badge variant="archive" title={recipe.docker_image}>
            docker
          </Badge>
        ) : (
          <span className="rg-muted">NA</span>
        ),
    },
    { key: 'inputs', header: 'Inputs', render: inputSummary },
  ];

  return (
    <section className="flex flex-col gap-4">
      <h2 className="text-xl font-semibold">Recipes</h2>
      <p className="rg-muted text-sm">
        Read only. Install recipes with{' '}
        <code className="rg-code rg-code--inline">refgenie recipe add</code>.
      </p>
      <DataTable
        caption="Recipes"
        columns={columns}
        rows={recipes.data?.results}
        rowKey={(recipe, index) => `${recipe.name}-${recipe.version}-${index}`}
        loading={recipes.isPending}
        error={recipes.error}
        onRetry={() => recipes.refetch()}
      />
    </section>
  );
}

export function AssetClassPanel() {
  const assetClasses = useAssetClasses({ pageSize: LIMIT });

  const columns: Array<Column<AssetClassPublic>> = [
    {
      key: 'name',
      header: 'Asset class',
      render: (cls) =>
        cls.id !== null && cls.id !== undefined ? (
          <Link className="rg-link font-medium" to={`/asset-classes/${cls.id}`}>
            {cls.name}
          </Link>
        ) : (
          cls.name
        ),
    },
    { key: 'version', header: 'Version', render: (cls) => cls.version },
    {
      key: 'description',
      header: 'Description',
      className: 'rg-table__cell--wrap',
      render: (cls) => cls.description ?? <span className="rg-muted">NA</span>,
    },
    {
      key: 'serving',
      header: 'Serving modes',
      render: (cls) => cls.serving_modes.join(', ') || <span className="rg-muted">NA</span>,
    },
  ];

  return (
    <section className="flex flex-col gap-4">
      <h2 className="text-xl font-semibold">Asset classes</h2>
      <p className="rg-muted text-sm">
        Read only. Install asset classes with{' '}
        <code className="rg-code rg-code--inline">refgenie asset_class add</code>.
      </p>
      <DataTable
        caption="Asset classes"
        columns={columns}
        rows={assetClasses.data?.results}
        rowKey={(cls, index) => `${cls.name}-${cls.version}-${index}`}
        loading={assetClasses.isPending}
        error={assetClasses.error}
        onRetry={() => assetClasses.refetch()}
      />
    </section>
  );
}
