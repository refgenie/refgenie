/**
 * `/species`: one row per organism, searchable by any name or taxon ID.
 * Built in the browser from the cached genome list (`useSpeciesIndex`), so it
 * works in both modes.
 */

import { useMemo } from 'react';
import { useSearchParamsState } from '../hooks/useSearchParamsState';
import { useCapability } from '../hooks/useCapability';
import { useRouteBase } from '../hooks/useRouteBase';
import { useSpeciesIndex } from '../hooks/queries/useSpecies';
import { MiniHero } from '../components/layout/MiniHero';
import { SpeciesTable } from '../components/species/SpeciesTable';
import { SearchBox } from '../components/common/SearchBox';
import { Pagination } from '../components/common/Pagination';
import { EmptyState } from '../components/common/states';
import { compareSpecies, speciesMatches } from '../utils/species';
import type { PageMeta } from '../types/pagination';

export function SpeciesPage() {
  const [state, actions] = useSearchParamsState();
  const index = useSpeciesIndex();

  // Call both hooks unconditionally; `&&` would make one conditional.
  const routeBase = useRouteBase();
  const canSeeTree = useCapability('archives');
  const showTreeLink = routeBase === '' && canSeeTree;

  const all = useMemo(() => [...(index.data?.values() ?? [])], [index.data]);

  const matched = useMemo(
    () => all.filter((row) => speciesMatches(row, state.q)).sort(compareSpecies),
    [all, state.q],
  );

  // Built locally so the shared <Pagination> works over an in-memory list.
  const pagination: PageMeta = {
    page: state.page,
    page_size: state.pageSize,
    total: matched.length,
  };
  const first = state.page * state.pageSize;
  const rows = index.data ? matched.slice(first, first + state.pageSize) : undefined;

  const named = all.filter((row) => row.key !== '');
  const genomeTotal = all.reduce((sum, row) => sum + row.genomes, 0);
  const withAssets = named.filter((row) => row.assets > 0).length;

  return (
    <div className="flex flex-col gap-6">
      <MiniHero
        title="Species"
        lede={
          <>
            Each species and how many genomes refgenie has for it.
            {index.data && (
              <>
                {' '}
                {named.length} species · {genomeTotal} genomes · {withAssets} species
                with built assets.
              </>
            )}
          </>
        }
      />

      <SearchBox
        value={state.q}
        onChange={actions.setQuery}
        placeholder="Search species by any name or taxon ID…"
        label="Search species"
      />

      <SpeciesTable
        species={rows}
        showTreeLink={showTreeLink}
        loading={index.isPending}
        error={index.error}
        onRetry={() => index.refetch()}
        empty={
          <EmptyState
            query={state.q || undefined}
            message="No genomes are registered here yet, so there are no species to list."
            onClearSearch={actions.clearSearch}
          />
        }
      />

      <Pagination pagination={pagination} onPageChange={actions.setPage} />
    </div>
  );
}
