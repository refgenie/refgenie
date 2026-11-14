/** `/collections/:id`: one branded collection, its genomes and its species tree. */

import { lazy, Suspense, useCallback, useId, useMemo } from 'react';
import { useParams, useSearchParams } from 'react-router-dom';
import { CollectionHero } from '../components/collections/CollectionHero';
import { CollectionStats } from '../components/collections/CollectionStats';
import { CollectionGenomeDetail } from '../components/collections/CollectionGenomeDetail';
import {
  COLLECTION_HIDDEN_COLUMNS,
  collectionColumns,
} from '../components/collections/collectionColumns';
import { useAccent } from '../components/collections/useAccent';
import { GenomeTable } from '../components/genomes/GenomeTable';
import { Pagination } from '../components/common/Pagination';
import { ScopedLink } from '../components/common/ScopedLink';
import { EmptyState, LoadingState } from '../components/common/states';
import { useDocumentTitle } from '../hooks/useDocumentTitle';
import { useCollection, useCollectionGenomes } from '../hooks/queries/useCollections';
import { collectionStats, filterCollectionGenomes } from '../utils/collections';
import { HAPLOTYPES } from '../utils/haplotype';
import type { Haplotype } from '../utils/haplotype';
import { DEFAULT_PAGE_SIZE } from '../types/pagination';
import type { GenomeResponse } from '../types/api';

const CollectionTree = lazy(() => import('../components/collections/CollectionTree'));

function asHaplotype(value: string | null): Haplotype | null {
  return value && (HAPLOTYPES as readonly string[]).includes(value)
    ? (value as Haplotype)
    : null;
}

function NotFound({ id }: { id: string | undefined }) {
  useDocumentTitle('Collection not found');
  return (
    <div className="flex flex-col gap-4">
      <EmptyState message={`No collection named “${id ?? ''}” on this server.`} />
      <ScopedLink className="rg-link" to="/collections">
        See all collections
      </ScopedLink>
    </div>
  );
}

export function CollectionPage() {
  const { id } = useParams<{ id: string }>();
  const collection = useCollection(id);
  if (!collection) return <NotFound id={id} />;
  // Keyed so a different collection starts with fresh state.
  return <CollectionView key={collection.id} id={collection.id} />;
}

function CollectionView({ id }: { id: string }) {
  const collection = useCollection(id)!;
  const query = useCollectionGenomes(collection);
  const accentRef = useAccent<HTMLDivElement>(collection.accent_color);
  const all = query.data;
  const [params, setParams] = useSearchParams();
  const searchId = useId();
  const haplotypeId = useId();
  const assetsId = useId();

  const q = params.get('q') ?? '';
  const haplotype = asHaplotype(params.get('haplotype'));
  const onlyWithAssets = params.get('assets') === 'built';
  const page = Math.max(0, Number(params.get('page') ?? 0) || 0);
  const selectedDigest = params.get('genome');

  const update = useCallback(
    (changes: Record<string, string | null>, resetPage = true) => {
      setParams(
        (previous) => {
          const next = new URLSearchParams(previous);
          for (const [key, value] of Object.entries(changes)) {
            if (value === null || value === '') next.delete(key);
            else next.set(key, value);
          }
          if (resetPage) next.delete('page');
          return next;
        },
        { replace: true },
      );
    },
    [setParams],
  );

  const stats = useMemo(() => (all ? collectionStats(all) : undefined), [all]);
  const filtered = useMemo(
    () =>
      all ? filterCollectionGenomes(all, { q, haplotype, onlyWithAssets }) : undefined,
    [all, q, haplotype, onlyWithAssets],
  );
  const pageRows = filtered?.slice(
    page * DEFAULT_PAGE_SIZE,
    (page + 1) * DEFAULT_PAGE_SIZE,
  );
  const selected = selectedDigest
    ? all?.find((g) => g.digest === selectedDigest)
    : undefined;

  const extraColumns = useMemo(
    () =>
      collectionColumns(collection, (genome: GenomeResponse) =>
        update({ genome: genome.digest }, false),
      ),
    [collection, update],
  );

  const showHaplotype = collection.columns.includes('haplotype');

  return (
    <div className="rg-collection flex flex-col gap-6" ref={accentRef}>
      <CollectionHero collection={collection} />
      <CollectionStats collection={collection} stats={stats} />

      <section className="flex flex-col gap-3" aria-labelledby={`${searchId}-heading`}>
        <h2 id={`${searchId}-heading`} className="text-xl font-semibold">
          Genomes
        </h2>
        <div className="rg-collection__filters flex flex-wrap items-center gap-3">
          <label className="sr-only" htmlFor={searchId}>
            Search this collection
          </label>
          <input
            id={searchId}
            className="rg-field__input rg-tree__filter"
            type="search"
            placeholder="Search by name, species or accession…"
            value={q}
            onChange={(event) => update({ q: event.target.value })}
          />
          {showHaplotype && (
            <>
              <label className="sr-only" htmlFor={haplotypeId}>
                Haplotype
              </label>
              <select
                id={haplotypeId}
                className="rg-field__input w-auto"
                value={haplotype ?? ''}
                onChange={(event) => update({ haplotype: event.target.value })}
              >
                <option value="">Any haplotype</option>
                {HAPLOTYPES.map((h) => (
                  <option key={h} value={h}>
                    {h}
                  </option>
                ))}
              </select>
            </>
          )}
          <label className="rg-check rg-check--inline" htmlFor={assetsId}>
            <input
              id={assetsId}
              type="checkbox"
              checked={onlyWithAssets}
              onChange={(event) =>
                update({ assets: event.target.checked ? 'built' : null })
              }
            />
            Only genomes with built assets
          </label>
        </div>

        <GenomeTable
          genomes={pageRows}
          query={q || undefined}
          extraColumns={extraColumns}
          hideColumns={COLLECTION_HIDDEN_COLUMNS}
          loading={query.isPending}
          error={query.error}
          onRetry={() => query.refetch()}
          empty={
            <EmptyState
              query={q || undefined}
              message="No genomes in this collection match these filters."
              onClearSearch={() => update({ q: null })}
            />
          }
        />
        <Pagination
          pagination={
            filtered
              ? { page, page_size: DEFAULT_PAGE_SIZE, total: filtered.length }
              : undefined
          }
          onPageChange={(next) => update({ page: String(next) }, false)}
        />
      </section>

      <section className="flex flex-col gap-3" aria-labelledby={`${searchId}-tree`}>
        <h2 id={`${searchId}-tree`} className="text-xl font-semibold">
          Species tree
        </h2>
        {all && (
          <Suspense fallback={<LoadingState label="Loading the tree" />}>
            <CollectionTree genomes={all} />
          </Suspense>
        )}
      </section>

      {selected && (
        <CollectionGenomeDetail
          collection={collection}
          genome={selected}
          onClose={() => update({ genome: null }, false)}
        />
      )}
    </div>
  );
}
