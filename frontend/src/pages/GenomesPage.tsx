import { useId } from 'react';
import { useSearchParams } from 'react-router-dom';
import { useSearchParamsState } from '../hooks/useSearchParamsState';
import { useResetPastEndPage } from '../hooks/useResetPastEndPage';
import { useGenomes } from '../hooks/queries/useGenomes';
import { useSummary } from '../hooks/queries/useServerInfo';
import { useCapability } from '../hooks/useCapability';
import { useRouteBase } from '../hooks/useRouteBase';
import { useBridgeDigests } from '../hooks/queries/useBridgeDigests';
import { useBridgeStore } from '../stores/bridgeStore';
import { MiniHero } from '../components/layout/MiniHero';
import { GenomeTable } from '../components/genomes/GenomeTable';
import { Badge } from '../components/common/Badge';
import { SearchBox } from '../components/common/SearchBox';
import { Pagination } from '../components/common/Pagination';
import { EmptyState } from '../components/common/states';
import { GENOME_SEARCH_FIELDS } from '../services/resources/genomes';
import type { Column } from '../components/common/DataTable';
import type { GenomeResponse } from '../types/api';

/** Most genomes have no assets, so the list hides them by default; `?assets=all` shows all. */
const WITH_ASSETS_PARAM = 'assets';
const SHOW_ALL = 'all';

export function GenomesPage() {
  const [state, actions] = useSearchParamsState();
  const [params, setParams] = useSearchParams();
  const canReadArchives = useCapability('archives');
  // The empty state suggests only what this instance can do.
  const canInitGenome = useCapability('genome_init');
  const canBuild = useCapability('build');
  const canPullRemote = useCapability('remote_browse');
  const noGenomesMessage =
    canInitGenome || canBuild
      ? `No genomes yet. Initialize one from a FASTA or build an asset from the Manage page${
          canPullRemote
            ? ', or pull one from a subscribed server on the Remote page'
            : ''
        }.`
      : 'No genomes are registered here yet.';
  const toggleId = useId();

  // Presence badges (is this genome also local?) are pointless under `/local`.
  const bridgeConnected = useBridgeStore((state) => state.status) === 'connected';
  const routeBase = useRouteBase();
  const bridgeDigests = useBridgeDigests(routeBase === '');
  const localColumn: Array<Column<GenomeResponse>> =
    bridgeConnected && routeBase === ''
      ? [
          {
            key: 'local',
            header: 'Local',
            render: (genome) =>
              bridgeDigests.data?.genomeDigests.has(genome.digest) ? (
                <Badge
                  variant="local"
                  title="You already have this genome on your local refgenie"
                >
                  local
                </Badge>
              ) : (
                <span className="rg-muted">not local</span>
              ),
          },
        ]
      : [];

  const onlyWithAssets = params.get(WITH_ASSETS_PARAM) !== SHOW_ALL;

  const query = useGenomes({
    has_assets: onlyWithAssets || undefined,
    q: state.q || undefined,
    searchFields: state.fields.length ? state.fields : undefined,
    operator: state.operator,
    page: state.page,
    pageSize: state.pageSize,
  });
  useResetPastEndPage(query.error, state.page, actions.setPage);

  const summary = useSummary({ enabled: canReadArchives });
  // Only the total is read.
  const withAssets = useGenomes(
    { has_assets: true, pageSize: 1 },
    { enabled: canReadArchives },
  );

  // When the filter empties a search, count the hidden matches for the empty state.
  const filteredOutSearch =
    onlyWithAssets && !!state.q && query.data?.pagination.total === 0;
  const unfilteredMatches = useGenomes(
    {
      q: state.q || undefined,
      searchFields: state.fields.length ? state.fields : undefined,
      operator: state.operator,
      pageSize: 1,
    },
    { enabled: filteredOutSearch },
  );
  const hiddenMatches = filteredOutSearch
    ? (unfilteredMatches.data?.pagination.total ?? 0)
    : 0;

  const setOnlyWithAssets = (on: boolean) => {
    setParams(
      (previous) => {
        const next = new URLSearchParams(previous);
        if (on) next.delete(WITH_ASSETS_PARAM);
        else next.set(WITH_ASSETS_PARAM, SHOW_ALL);
        // Page 4 of the filtered list is not page 4 of the full one.
        next.delete('page');
        return next;
      },
      { replace: true },
    );
  };

  return (
    <div className="flex flex-col gap-6">
      <MiniHero
        title="Genomes"
        lede={
          <>
            A genome here is one reference assembly, identified by a digest computed
            from its sequences rather than by whatever name a provider gave it. Each row
            lists the names that genome answers to and how many assets have been built
            for it.
          </>
        }
      />

      {summary.data && (
        <dl className="grid grid-cols-1 sm:grid-cols-3 gap-4">
          <div className="rg-card p-4">
            <dt className="rg-kv__term">Genomes</dt>
            <dd className="text-2xl font-semibold">
              {withAssets.data ? withAssets.data.pagination.total : '—'}
              <span className="rg-muted text-sm font-normal">
                {' '}
                with assets · {summary.data.genomes} registered
              </span>
            </dd>
          </div>
          {(
            [
              ['Asset groups', summary.data.asset_groups],
              ['Assets', summary.data.assets],
            ] as const
          ).map(([term, value]) => (
            <div className="rg-card p-4" key={term}>
              <dt className="rg-kv__term">{term}</dt>
              <dd className="text-2xl font-semibold">{value}</dd>
            </div>
          ))}
        </dl>
      )}

      <SearchBox
        value={state.q}
        onChange={actions.setQuery}
        fields={GENOME_SEARCH_FIELDS}
        selectedFields={state.fields}
        onFieldsChange={actions.setFields}
        operator={state.operator}
        onOperatorChange={actions.setOperator}
        placeholder="Search genomes…"
        label="Search genomes"
      />

      <label className="rg-check rg-check--inline" htmlFor={toggleId}>
        <input
          id={toggleId}
          type="checkbox"
          checked={onlyWithAssets}
          onChange={(event) => setOnlyWithAssets(event.target.checked)}
        />
        Only genomes with assets
      </label>

      <GenomeTable
        genomes={query.data?.results}
        extraColumns={localColumn}
        query={state.q || undefined}
        loading={query.isPending}
        error={query.error}
        onRetry={() => query.refetch()}
        empty={
          hiddenMatches > 0 ? (
            <div className="rg-state rg-state--empty">
              <p className="mb-4">
                {hiddenMatches}{' '}
                {hiddenMatches === 1 ? 'genome matches' : 'genomes match'} “{state.q}”,
                but {hiddenMatches === 1 ? 'it has no' : 'none have'} assets built here.
              </p>
              <button
                type="button"
                className="rg-btn rg-btn--sm"
                onClick={() => setOnlyWithAssets(false)}
              >
                Show all matches
              </button>
            </div>
          ) : (
            <EmptyState
              query={state.q || undefined}
              message={
                onlyWithAssets
                  ? 'No genome here has any assets built. Untick the filter to see every registered genome.'
                  : noGenomesMessage
              }
              onClearSearch={actions.clearSearch}
            />
          )
        }
      />

      <Pagination pagination={query.data?.pagination} onPageChange={actions.setPage} />
    </div>
  );
}
