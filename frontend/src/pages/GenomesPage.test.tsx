import { describe, expect, it } from 'vitest';
import { screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { HttpResponse, http } from 'msw';
import { GenomesPage } from './GenomesPage';
import { server } from '../test/server';
import { API } from '../test/handlers';
import { emptyPage } from '../test/fixtures';
import { renderWithProviders } from '../test/renderWithProviders';

describe('GenomesPage', () => {
  it('renders the parity fields once data arrives', async () => {
    renderWithProviders(<GenomesPage />);
    expect(await screen.findByRole('link', { name: 'hg38' })).toBeInTheDocument();
    expect(screen.getByText('Human GRCh38 reference assembly')).toBeInTheDocument();
    expect(screen.getAllByText('Homo sapiens').length).toBeGreaterThan(0);
    expect(screen.getByText('NCBI · GCA_000001405.15')).toBeInTheDocument();
  });

  it('sends the search term to the server rather than filtering in memory', async () => {
    const seen: string[] = [];
    server.use(
      http.get(`${API}/genomes`, ({ request }) => {
        seen.push(new URL(request.url).search);
        return HttpResponse.json(emptyPage);
      }),
    );

    renderWithProviders(<GenomesPage />);
    await userEvent.type(screen.getByLabelText('Search genomes'), 'hg38');

    await waitFor(() => expect(seen.some((s) => s.includes('q=hg38'))).toBe(true));
  });

  it('renders the no-data empty state', async () => {
    server.use(http.get(`${API}/genomes`, () => HttpResponse.json(emptyPage)));
    renderWithProviders(<GenomesPage />);
    expect(await screen.findByText(/No genomes yet/)).toBeInTheDocument();
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
