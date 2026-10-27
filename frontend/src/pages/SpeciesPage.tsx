/**
 * `/species` — one row per organism, searchable by scientific name, common
 * name, or taxon ID in a single box.
 *
 * Derived entirely from the cached genome list (`useSpeciesIndex`), so this
 * page works on a local dash as well as on the server and needs no capability
 * gate. Filtering, sorting and paging happen over that in-memory slice, the
 * same way `GenomesPage` handles its assets filter.
 *
 * There is no "only species with assets" toggle. A species table is a
 * directory: hiding rows would make the header count a lie, and the assets-desc
 * default already floats the useful ones to the top.
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
import type { PaginationMeta } from '../types/pagination';

export function SpeciesPage() {
  const [state, actions] = useSearchParamsState();
  const index = useSpeciesIndex();

  // Both hooks run unconditionally and are combined afterwards. Writing this as
  // `useRouteBase() === '' && useCapability('archives')` would short-circuit
  // past a hook call — a conditional hook, and a lint error.
  const routeBase = useRouteBase();
  const canSeeTree = useCapability('archives');
  const showTreeLink = routeBase === '' && canSeeTree;

  const all = useMemo(() => [...(index.data?.values() ?? [])], [index.data]);

  const matched = useMemo(
    () => all.filter((row) => speciesMatches(row, state.q)).sort(compareSpecies),
    [all, state.q],
  );

  // Synthesized so the shared <Pagination> works unchanged over a local slice.
  const pagination: PaginationMeta = {
    offset: state.offset,
    limit: state.limit,
    total: matched.length,
  };
  const rows = index.data ? matched.slice(state.offset, state.offset + state.limit) : undefined;

  const named = all.filter((row) => row.key !== '');
  const withAssets = named.filter((row) => row.assets > 0).length;

  return (
    <div className="flex flex-col gap-6">
      <MiniHero
        title="Species"
        lede={
          <>
            One row per organism, gathering every genome refgenie holds for it. A species can
            have several assemblies — different builds, different providers — and they all fold
            into its one row here.
            {index.data && (
              <>
                {' '}
                {named.length} species · {withAssets} with assets built.
              </>
            )}
          </>
        }
      />

      {/*
        No field selector and no operator selector: the rows are derived in the
        browser, not fetched, so those controls would be theatre over a
        two-field object. One box matches both name kinds and the taxon ID.
      */}
      <SearchBox
        value={state.q}
        onChange={actions.setQuery}
        placeholder="Search species — scientific name, common name, or taxon ID…"
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

      <Pagination pagination={pagination} onOffsetChange={actions.setOffset} />
    </div>
  );
}
