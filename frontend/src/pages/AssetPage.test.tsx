import { describe, expect, it } from 'vitest';
import { screen } from '@testing-library/react';
import { HttpResponse, http } from 'msw';
import { Route, Routes } from 'react-router-dom';
import { AssetPage } from './AssetPage';
import { server } from '../test/server';
import { API } from '../test/handlers';
import { renderWithProviders, serverConfig } from '../test/renderWithProviders';
import type { UiConfig } from '../types/ui';

const DIGEST = '1111111111111111111111111111111111111111111111111111111111111111';

function renderAsset(config?: UiConfig) {
  return renderWithProviders(
    <Routes>
      <Route path="/assets/:digest" element={<AssetPage />} />
    </Routes>,
    { route: `/assets/${DIGEST}`, config },
  );
}

describe('AssetPage', () => {
  it('renders the asset parity fields, including seek keys', async () => {
    renderAsset();
    expect(await screen.findByRole('heading', { name: 'fasta:default', level: 1 })).toBeInTheDocument();
    expect(screen.getByText('Seek keys')).toBeInTheDocument();
    // The seek key row: name, type badge, value.
    expect(screen.getAllByText('hg38.fa').length).toBeGreaterThan(0);
    expect(screen.getByText(/The FASTA file/)).toBeInTheDocument();
  });

  it('lists parent assets from the expanded relationships endpoint', async () => {
    renderAsset();
    expect(await screen.findByText('Parents')).toBeInTheDocument();
    expect(screen.getByText('No child assets.')).toBeInTheDocument();
  });

  it('links files for download only when the downloads capability is on', async () => {
    renderAsset(serverConfig);
    const link = await screen.findByRole('link', { name: 'hg38.fa.fai' });
    expect(link).toHaveAttribute('href', `/v4/assets/${DIGEST}/files/hg38.fa.fai`);
  });

  it('renders ErrorState on a 500', async () => {
    server.use(
      http.get(`${API}/assets/:digest`, () =>
        HttpResponse.json({ detail: 'kaboom' }, { status: 500 }),
      ),
    );
    renderAsset();
    expect(await screen.findByRole('alert')).toHaveTextContent('kaboom');
  });
});
