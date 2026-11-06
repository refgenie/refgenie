import { describe, expect, it } from 'vitest';
import { screen, waitFor, within } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { HttpResponse, http } from 'msw';
import { GenomesPage } from './GenomesPage';
import { server } from '../test/server';
import { API } from '../test/handlers';
import { emptyPage, genomesFixture } from '../test/fixtures';
import { renderWithProviders } from '../test/renderWithProviders';
import type { GenomeResponse } from '../types/api';

/** A registered genome with nothing built — 675 of the public server's 701. */
const BARE: GenomeResponse = {
  digest: '0000000000000000000000000000beef',
  aliases: ['bare-genome'],
  description: null,
  asset_count: 0,
  species_name: null,
  common_name: null,
  synonyms: null,
  taxon_id: null,
  assembly_source: null,
  assembly_accession: null,
  species: null,
};

/** Answers like the server, including its `has_assets` filter. */
function seed(items: GenomeResponse[]) {
  server.use(
    http.get(`${API}/genomes`, ({ request }) => {
      const withAssets = new URL(request.url).searchParams.get('has_assets') === 'true';
      const results = withAssets
        ? items.filter((genome) => genome.asset_count > 0)
        : items;
      return HttpResponse.json({
        results,
        pagination: { page: 0, page_size: 1000, total: results.length },
      });
    }),
  );
}

describe('GenomesPage', () => {
  it('renders the parity fields once data arrives', async () => {
    renderWithProviders(<GenomesPage />);
    expect(await screen.findByRole('link', { name: 'hg38' })).toBeInTheDocument();
    expect(screen.getByText('Human GRCh38 reference assembly')).toBeInTheDocument();
    expect(screen.getAllByText('Homo sapiens').length).toBeGreaterThan(0);
    expect(screen.getByText('NCBI · GCA_000001405.15')).toBeInTheDocument();
  });

  it('shows one common name beside the scientific name', async () => {
    renderWithProviders(<GenomesPage />);
    const row = await screen.findByRole('row', { name: /hg38/ });
    expect(within(row).getByText('human')).toHaveClass('rg-muted');
    expect(row).toHaveTextContent('Homo sapiens human');
  });

  it('hides genomes with no assets by default', async () => {
    seed([...genomesFixture.results, BARE]);
    renderWithProviders(<GenomesPage />);
    await screen.findByRole('link', { name: 'hg38' });
    expect(screen.queryByRole('link', { name: 'bare-genome' })).not.toBeInTheDocument();
  });

  it('shows them once the filter is unticked, and records that in the URL', async () => {
    seed([...genomesFixture.results, BARE]);
    renderWithProviders(<GenomesPage />);
    await screen.findByRole('link', { name: 'hg38' });

    await userEvent.click(screen.getByLabelText('Only genomes with assets'));
    expect(
      await screen.findByRole('link', { name: 'bare-genome' }),
    ).toBeInTheDocument();
  });

  it('keeps a deep link with ?assets=all unfiltered', async () => {
    seed([...genomesFixture.results, BARE]);
    renderWithProviders(<GenomesPage />, { route: '/genomes?assets=all' });
    expect(
      await screen.findByRole('link', { name: 'bare-genome' }),
    ).toBeInTheDocument();
  });

  it('sends the search term and the assets filter to the server', async () => {
    const seen: URLSearchParams[] = [];
    server.use(
      http.get(`${API}/genomes`, ({ request }) => {
        seen.push(new URL(request.url).searchParams);
        return HttpResponse.json(emptyPage);
      }),
    );

    renderWithProviders(<GenomesPage />);
    await userEvent.type(screen.getByLabelText('Search genomes'), 'rCRSd');

    await waitFor(() =>
      expect(
        seen.some((s) => s.get('q') === 'rCRSd' && s.get('has_assets') === 'true'),
      ).toBe(true),
    );
  });

  it('sends the search term to the server when the assets filter is off', async () => {
    const seen: string[] = [];
    server.use(
      http.get(`${API}/genomes`, ({ request }) => {
        seen.push(new URL(request.url).search);
        return HttpResponse.json(emptyPage);
      }),
    );

    renderWithProviders(<GenomesPage />, { route: '/genomes?assets=all' });
    await userEvent.type(screen.getByLabelText('Search genomes'), 'hg38');

    await waitFor(() => expect(seen.some((s) => s.includes('q=hg38'))).toBe(true));
  });

  it('searches the species field, which covers every species name', async () => {
    const seen: URLSearchParams[] = [];
    server.use(
      http.get(`${API}/genomes`, ({ request }) => {
        seen.push(new URL(request.url).searchParams);
        return HttpResponse.json(emptyPage);
      }),
    );

    renderWithProviders(<GenomesPage />);
    await userEvent.selectOptions(screen.getByLabelText('Field to search'), 'species');
    await userEvent.type(screen.getByLabelText('Search genomes'), 'brown rat');

    await waitFor(() =>
      expect(
        seen.some(
          (s) => s.get('q') === 'brown rat' && s.get('search_fields') === 'species',
        ),
      ).toBe(true),
    );
  });

  it('renders the no-data empty state', async () => {
    server.use(http.get(`${API}/genomes`, () => HttpResponse.json(emptyPage)));
    renderWithProviders(<GenomesPage />, { route: '/genomes?assets=all' });
    expect(await screen.findByText(/No genomes yet/)).toBeInTheDocument();
  });

  it('names the filter when it is what emptied the list', async () => {
    seed([BARE]);
    renderWithProviders(<GenomesPage />);
    expect(await screen.findByText(/Untick the filter/)).toBeInTheDocument();
  });

  it('counts the matches the filter hid and offers to show them', async () => {
    seed([BARE]);
    renderWithProviders(<GenomesPage />, { route: '/genomes?q=bare' });
    expect(
      await screen.findByText(/1 genome matches “bare”, but it has no assets/),
    ).toBeInTheDocument();
    await userEvent.click(screen.getByRole('button', { name: 'Show all matches' }));
    expect(
      await screen.findByRole('link', { name: 'bare-genome' }),
    ).toBeInTheDocument();
  });
});
