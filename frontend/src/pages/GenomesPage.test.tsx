import { describe, expect, it } from 'vitest';
import { screen, waitFor } from '@testing-library/react';
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
  taxon_id: null,
  assembly_source: null,
  assembly_accession: null,
};

/** Both endpoints the page can be driven from answer with the same rows. */
function seed(items: GenomeResponse[]) {
  server.use(
    http.get(`${API}/genomes`, () =>
      HttpResponse.json({
        items,
        pagination: { offset: 0, limit: 1000, total: items.length },
      }),
    ),
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

  it('hides genomes with no assets by default', async () => {
    seed([...genomesFixture.items, BARE]);
    renderWithProviders(<GenomesPage />);
    await screen.findByRole('link', { name: 'hg38' });
    expect(screen.queryByRole('link', { name: 'bare-genome' })).not.toBeInTheDocument();
  });

  it('shows them once the filter is unticked, and records that in the URL', async () => {
    seed([...genomesFixture.items, BARE]);
    renderWithProviders(<GenomesPage />);
    await screen.findByRole('link', { name: 'hg38' });

    await userEvent.click(screen.getByLabelText('Only genomes with assets'));
    expect(await screen.findByRole('link', { name: 'bare-genome' })).toBeInTheDocument();
  });

  it('keeps a deep link with ?assets=all unfiltered', async () => {
    seed([...genomesFixture.items, BARE]);
    renderWithProviders(<GenomesPage />, { route: '/genomes?assets=all' });
    expect(await screen.findByRole('link', { name: 'bare-genome' })).toBeInTheDocument();
  });

  it('filters in memory while the assets filter is on', async () => {
    seed([...genomesFixture.items, BARE]);
    renderWithProviders(<GenomesPage />);
    await screen.findByRole('link', { name: 'hg38' });

    await userEvent.type(screen.getByLabelText('Search genomes'), 'rCRSd');
    await waitFor(() =>
      expect(screen.queryByRole('link', { name: 'hg38' })).not.toBeInTheDocument(),
    );
    expect(screen.getByRole('link', { name: 'rCRSd' })).toBeInTheDocument();
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

  it('titles the tab with the page and the service, not a bare "refgenie"', async () => {
    renderWithProviders(<GenomesPage />);
    await screen.findByRole('link', { name: 'hg38' });
    expect(document.title).toBe('Genomes · refgenie local dashboard');
  });

  it('renders ErrorState on a 500', async () => {
    server.use(
      http.get(`${API}/genomes`, () =>
        HttpResponse.json({ detail: 'boom' }, { status: 500 }),
      ),
    );
    renderWithProviders(<GenomesPage />);
    expect(await screen.findByRole('alert')).toHaveTextContent('boom');
  });
});
