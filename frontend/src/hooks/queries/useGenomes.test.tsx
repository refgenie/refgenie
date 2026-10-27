/**
 * The sharpest cache property in the app: the whole-backend genome list is ONE
 * entry, read through two different `select`s, and it costs ONE request per
 * session no matter how many pages read it.
 *
 * `/genomes` and `/species` are the two consumers; before this guarantee the
 * species page re-fetched 701 rows and rebuilt a Map that the genome page had
 * already paid for.
 */

import { describe, expect, it } from 'vitest';
import { HttpResponse, http } from 'msw';
import { render, screen } from '@testing-library/react';
import { useGenomeIndex } from './useGenomes';
import { useSpeciesIndex } from './useSpecies';
import { server } from '../../test/server';
import { genomesFixture } from '../../test/fixtures';
import { createWrapper } from '../../test/renderWithProviders';
import type { GenomeIndex } from '../../utils/genomes';
import type { SpeciesIndex } from '../../utils/species';

const genomeMaps: (GenomeIndex | undefined)[] = [];
const speciesMaps: (SpeciesIndex | undefined)[] = [];

function GenomeProbe() {
  const index = useGenomeIndex();
  genomeMaps.push(index.data);
  return <span data-testid="genomes">{index.data ? `genomes:${index.data.size}` : 'pending'}</span>;
}

function SpeciesProbe() {
  const index = useSpeciesIndex();
  speciesMaps.push(index.data);
  return <span data-testid="species">{index.data ? `species:${index.data.size}` : 'pending'}</span>;
}

/** `Array.prototype.at` is above this project's lib target. */
function last<T>(values: T[]): T | undefined {
  return values[values.length - 1];
}

/** Counts every hit on the list endpoint `fetchAllGenomes` pages through. */
function countGenomeRequests() {
  const counter = { hits: 0 };
  server.use(
    http.get('/v4/genomes', () => {
      counter.hits += 1;
      return HttpResponse.json(genomesFixture);
    }),
  );
  return counter;
}

describe('the genome index', () => {
  it('costs one request for both readers, and builds each Map once', async () => {
    genomeMaps.length = 0;
    speciesMaps.length = 0;
    const counter = countGenomeRequests();
    const Wrapper = createWrapper();

    // /genomes
    const genomesView = render(<GenomeProbe />, { wrapper: Wrapper });
    expect(await screen.findByText('genomes:2')).toBeInTheDocument();
    const firstGenomeMap = last(genomeMaps);
    genomesView.unmount();

    // …then /species, on the same cache. No second request.
    const speciesView = render(<SpeciesProbe />, { wrapper: Wrapper });
    expect(await screen.findByText('species:1')).toBeInTheDocument();
    const firstSpeciesMap = last(speciesMaps);
    speciesView.unmount();

    // …then back to /genomes.
    render(<GenomeProbe />, { wrapper: Wrapper });
    expect(await screen.findByText('genomes:2')).toBeInTheDocument();

    expect(counter.hits).toBe(1);
    // Same instance on the way back: the transform was memoized, not redone.
    expect(last(genomeMaps)).toBe(firstGenomeMap);
    expect(firstSpeciesMap).not.toBe(firstGenomeMap);
  });

  it('shares one request between two readers mounted together', async () => {
    genomeMaps.length = 0;
    speciesMaps.length = 0;
    const counter = countGenomeRequests();

    render(
      <>
        <GenomeProbe />
        <SpeciesProbe />
      </>,
      { wrapper: createWrapper() },
    );

    expect(await screen.findByText('genomes:2')).toBeInTheDocument();
    expect(await screen.findByText('species:1')).toBeInTheDocument();
    expect(counter.hits).toBe(1);
  });

  it('gives every reader of one select the same Map instance', async () => {
    genomeMaps.length = 0;
    countGenomeRequests();

    render(
      <>
        <GenomeProbe />
        <GenomeProbe />
      </>,
      { wrapper: createWrapper() },
    );

    expect((await screen.findAllByText('genomes:2')).length).toBe(2);
    const loaded = genomeMaps.filter((map): map is GenomeIndex => map !== undefined);
    expect(new Set(loaded).size).toBe(1);
  });
});
