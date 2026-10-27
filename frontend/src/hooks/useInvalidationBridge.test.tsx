/**
 * The bus -> bridge -> cache chain.
 *
 * `jobStore` and the action buttons publish DOMAIN keys ('genomes'), never
 * cache keys; this is the one place that turns them into a prefix match, and a
 * mounted list has to refresh when it fires.
 */

import { describe, expect, it } from 'vitest';
import { act, render, screen } from '@testing-library/react';
import { HttpResponse, http } from 'msw';
import { QUERY_PREFIXES, useInvalidationBridge } from './useInvalidationBridge';
import { useGenomeIndex, useGenomes } from './queries/useGenomes';
import { useSpeciesIndex } from './queries/useSpecies';
import { useRecipes } from './queries/useRecipes';
import { server } from '../test/server';
import { emptyPage, genomesFixture } from '../test/fixtures';
import { createWrapper } from '../test/renderWithProviders';
import { invalidate } from '../services/invalidation';
import { qk } from '../services/queryKeys';

function Screen() {
  useInvalidationBridge();
  const genomes = useGenomes({ limit: 20 });
  const recipes = useRecipes({ limit: 20 });
  return (
    <>
      <span data-testid="genomes">{genomes.data?.items.length ?? 'pending'}</span>
      <span data-testid="recipes">{recipes.data?.items.length ?? 'pending'}</span>
    </>
  );
}

/** `/genomes` (its default view) and `/species` render from this one entry. */
function IndexScreen() {
  useInvalidationBridge();
  const genomes = useGenomeIndex();
  const species = useSpeciesIndex();
  return (
    <>
      <span data-testid="index">{genomes.data?.size ?? 'pending'}</span>
      <span data-testid="species">{species.data?.size ?? 'pending'}</span>
    </>
  );
}

describe('useInvalidationBridge', () => {
  it('refreshes the lists a domain key covers, and only those', async () => {
    let genomeHits = 0;
    let recipeHits = 0;
    server.use(
      http.get('/v4/genomes', () => {
        genomeHits += 1;
        return HttpResponse.json(genomesFixture);
      }),
      http.get('/v4/recipes', () => {
        recipeHits += 1;
        return HttpResponse.json(emptyPage);
      }),
    );

    render(<Screen />, { wrapper: createWrapper() });
    expect(await screen.findByText('2')).toBeInTheDocument();
    expect(genomeHits).toBe(1);
    expect(recipeHits).toBe(1);

    await act(async () => {
      invalidate(['genomes']);
      await new Promise((resolve) => setTimeout(resolve, 0));
    });

    expect(genomeHits).toBe(2);
    expect(recipeHits).toBe(1);
  });

  it('refreshes the whole-backend genome index, which is all /genomes and /species read', async () => {
    // The index is fetched with the bulk page size, so this handler answers
    // both `useGenomeIndex` and `useSpeciesIndex` from the same entry.
    let hits = 0;
    let rows = genomesFixture.items;
    server.use(
      http.get('/v4/genomes', () => {
        hits += 1;
        return HttpResponse.json({
          items: rows,
          pagination: { offset: 0, limit: 1000, total: rows.length },
        });
      }),
    );

    render(<IndexScreen />, { wrapper: createWrapper() });
    expect(await screen.findByTestId('index')).toHaveTextContent('2');
    expect(hits).toBe(1);

    // A pull finished and added a genome. `jobStore` publishes 'genomes'.
    rows = [...genomesFixture.items, { ...genomesFixture.items[0], digest: 'pulled-digest' }];

    await act(async () => {
      invalidate(['genomes']);
      await new Promise((resolve) => setTimeout(resolve, 0));
    });

    expect(hits).toBe(2);
    expect(await screen.findByTestId('index')).toHaveTextContent('3');
  });

  it('names every cache-key family exactly once, and none that no longer exists', () => {
    /**
     * Families deliberately off the domain bus, each for its own reason:
     *  - `ui-config` is resolved at startup and never changes under a running app.
     *  - `archives` exists only in server mode, where nothing publishes here.
     *  - `job-log` is paged by line offset; a blanket invalidation would refetch
     *    every open log, and `JobDetailModal` refetches on open anyway.
     *  - `bridge-digests` describes a DIFFERENT machine. Invalidating it from
     *    this backend's actions is the bug `useBridgeDigests` documents.
     */
    const OFF_THE_BUS = new Set(['ui-config', 'archives', 'job-log', 'bridge-digests']);

    const families = new Set(
      Object.values(qk).map(
        (make) => (make as (...args: never[]) => readonly unknown[])()[0] as string,
      ),
    );
    const mapped = new Set(Object.values(QUERY_PREFIXES).flat());

    // Nothing invalidatable may be forgotten...
    expect([...families].filter((f) => !OFF_THE_BUS.has(f) && !mapped.has(f))).toEqual([]);
    // ...and no prefix may outlive the key it was written for.
    expect([...mapped].filter((f) => !families.has(f))).toEqual([]);
  });
});
