import { useId, useMemo } from 'react';
import { useSearchParams } from 'react-router-dom';
import { useSearchParamsState } from '../hooks/useSearchParamsState';
import { useGenomeIndex, useGenomes } from '../hooks/queries/useGenomes';
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
import { genomeMatches, preferredAlias } from '../utils/genomes';
import type { Column } from '../components/common/DataTable';
import type { GenomeResponse } from '../types/api';
import type { PaginationMeta } from '../types/pagination';

/**
 * 675 of the public server's 701 genomes have nothing built. Paging through
 * fifteen screens of zeroes is not browsing, so the list hides them by default.
 *
 * `?assets=all` turns the filter off. The param is read straight off the URL
 * rather than through `useSearchParamsState`: five pages share that hook and
 * none of the others has this axis.
 */
const WITH_ASSETS_PARAM = 'assets';
const SHOW_ALL = 'all';

export function GenomesPage() {
  const [state, actions] = useSearchParamsState();
  const [params, setParams] = useSearchParams();
  const canReadArchives = useCapability('archives');
  // The empty state points at what THIS instance can do, not at the CLI: a
  // local dash initializes and builds from Manage and pulls from Remote; the
  // public site can do none of those.
  const canInitGenome = useCapability('genome_init');
  const canBuild = useCapability('build');
  const canPullRemote = useCapability('remote_browse');
  const noGenomesMessage = canInitGenome || canBuild
    ? `No genomes yet. Initialize one from a FASTA or build an asset from the Manage page${
        canPullRemote ? ', or pull one from a subscribed server on the Remote page' : ''
      }.`
    : 'No genomes are registered here yet.';
  const toggleId = useId();

  // Presence badges only make sense on the remote branch: digests are
  // content-derived, so remote-digest ∈ local-index is an exact match. Inside
  // `/local` every row is local by definition, so the column is dropped.
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

  // No endpoint has a `has_assets` parameter, or any sort parameter, so the
  // filtered view is driven from the cached full-list index and filtered,
  // sorted and paginated here. The unfiltered view stays server-paginated.
  const index = useGenomeIndex();
  const query = useGenomes(
    {
      q: state.q || undefined,
      searchFields: state.fields.length ? state.fields : undefined,
      operator: state.operator,
      offset: state.offset,
      limit: state.limit,
    },
    { enabled: !onlyWithAssets },
  );

  // Fills the <h3>Summary</h3> heading that was empty in the Jinja page.
  const summary = useSummary({ enabled: canReadArchives });

  const withAssets = useMemo(
    () => [...(index.data?.values() ?? [])].filter((genome) => genome.asset_count > 0),
    [index.data],
  );

  const matched = useMemo(() => {
    const label = (genome: GenomeResponse) =>
      preferredAlias(genome.aliases, state.q) ?? genome.digest;
    return withAssets
      .filter((genome) => genomeMatches(genome, state.q, state.fields, state.operator))
      .sort((a, b) => b.asset_count - a.asset_count || label(a).localeCompare(label(b)));
  }, [withAssets, state.q, state.fields, state.operator]);

  // Synthesized so the shared <Pagination> works unchanged over a local slice.
  const clientPagination: PaginationMeta = {
    offset: state.offset,
    limit: state.limit,
    total: matched.length,
  };

  const rows = onlyWithAssets
    ? matched.slice(state.offset, state.offset + state.limit)
    : query.data?.items;
  const pagination = onlyWithAssets ? clientPagination : query.data?.pagination;
  const source = onlyWithAssets ? index : query;

  const setOnlyWithAssets = (on: boolean) => {
    setParams(
      (previous) => {
        const next = new URLSearchParams(previous);
        if (on) next.delete(WITH_ASSETS_PARAM);
        else next.set(WITH_ASSETS_PARAM, SHOW_ALL);
        // Page 4 of the filtered list is not page 4 of the full one.
        next.delete('offset');
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
            A genome here is one reference assembly, identified by a digest computed from its
            sequences rather than by whatever name a provider gave it. Each row lists the names
            that genome answers to and how many assets have been built for it.
          </>
        }
      />

      {summary.data && (
        <dl className="grid grid-cols-1 sm:grid-cols-3 gap-4">
          <div className="rg-card p-4">
            <dt className="rg-kv__term">Genomes</dt>
            {/* The registered count is not the number of things anyone can
                download. Both numbers, or the number is a lie. */}
            <dd className="text-2xl font-semibold">
              {index.data ? withAssets.length : '—'}
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
        genomes={rows}
        extraColumns={localColumn}
        query={state.q || undefined}
        loading={source.isPending}
        error={source.error}
        onRetry={() => source.refetch()}
        empty={
          <EmptyState
            query={state.q || undefined}
            message={
              onlyWithAssets
                ? 'No genome here has any assets built. Untick the filter to see every registered genome.'
                : noGenomesMessage
            }
            onClearSearch={actions.clearSearch}
          />
        }
      />

      <Pagination pagination={pagination} onOffsetChange={actions.setOffset} />
    </div>
  );
}
