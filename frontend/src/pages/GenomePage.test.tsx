import { describe, expect, it } from 'vitest';
import { screen } from '@testing-library/react';
import { HttpResponse, http } from 'msw';
import { Route, Routes } from 'react-router-dom';
import { GenomePage } from './GenomePage';
import { server } from '../test/server';
import { API } from '../test/handlers';
import { renderWithProviders } from '../test/renderWithProviders';

const DIGEST = '5cb5b9d2d1e2b5e5e2f0a3c4d5e6f708';

function renderGenome(route = `/genomes/${DIGEST}`) {
  return renderWithProviders(
    <Routes>
      <Route path="/genomes/:digest" element={<GenomePage />} />
    </Routes>,
    { route },
  );
}

describe('GenomePage', () => {
  it('renders the genome detail parity fields', async () => {
    renderGenome();
    expect(await screen.findByRole('heading', { name: 'hg38', level: 1 })).toBeInTheDocument();
    expect(screen.getByText('chromosome')).toBeInTheDocument();
    // Taxon and accession are rendered as external links. The name carries the
    // sr-only "(opens in a new tab)" suffix that ExternalLink appends, because
    // the arrow icon itself is aria-hidden.
    expect(screen.getByRole('link', { name: /^9606\b/ })).toHaveAttribute(
      'href',
      'https://identifiers.org/taxonomy:9606',
    );
    expect(screen.getByRole('link', { name: /^GCA_000001405\.15\b/ })).toHaveAttribute(
      'href',
      'https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000001405.15/',
    );
  });

  it('shows the FHR panel', async () => {
    renderGenome();
    expect(await screen.findByText('FHR metadata')).toBeInTheDocument();
  });

  it('lists file-mode assets that have no archive record', async () => {
    renderGenome();
    // fasta_index is serving_modes: ["file"] and has no tarball, but the table
    // is driven from /assets so it still gets a row.
    expect(await screen.findByText('FASTA index')).toBeInTheDocument();
  });

  it('renders 404 copy naming the digest', async () => {
    server.use(
      http.get(`${API}/genomes/:digest`, () =>
        HttpResponse.json({ detail: 'Genome not found' }, { status: 404 }),
      ),
    );
    renderGenome();
    const alert = await screen.findByRole('alert');
    expect(alert).toHaveTextContent('Not found');
    expect(alert).toHaveTextContent(DIGEST);
  });

  // The head is a sibling of the loading/error branches, not something they
  // return past. A missing genome must still say where you are and offer a way
  // back, and must never title the tab with the "Genome" placeholder.
  it('keeps the head, and the way out, when the genome is missing', async () => {
    server.use(
      http.get(`${API}/genomes/:digest`, () =>
        HttpResponse.json({ detail: 'Genome not found' }, { status: 404 }),
      ),
    );
    renderGenome();
    await screen.findByRole('alert');
    expect(screen.getByRole('heading', { level: 1 })).toHaveTextContent('Genome');
    expect(screen.getByRole('link', { name: 'Genomes' })).toHaveAttribute('href', '/genomes');
    expect(document.title).toBe('refgenie local dashboard');
  });
});
