import { describe, expect, it, vi } from 'vitest';
import { screen, waitFor, within } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { HttpResponse, http } from 'msw';
import { Route, Routes } from 'react-router-dom';
import { AssetPage } from './AssetPage';
import { server } from '../test/server';
import { API } from '../test/handlers';
import { assetsFixture } from '../test/fixtures';
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

  it("names the asset's genome in the breadcrumb", async () => {
    renderAsset();
    const crumbs = await screen.findByRole('navigation', { name: /breadcrumb/i });
    expect(await within(crumbs).findByRole('link', { name: 'hg38' })).toHaveAttribute(
      'href',
      '/genomes/5cb5b9d2d1e2b5e5e2f0a3c4d5e6f708',
    );
  });

  it('lists parent assets from the expanded relationships endpoint', async () => {
    renderAsset();
    expect(await screen.findByText('Parents')).toBeInTheDocument();
    expect(screen.getByText('No child assets.')).toBeInTheDocument();
  });

  it('links files for download only when the downloads capability is on', async () => {
    renderAsset(serverConfig);
    // DownloadLink prefixes an sr-only "Download " so the accessible name puts
    // the verb first; the arrow icon itself is aria-hidden.
    const link = await screen.findByRole('link', { name: 'Download hg38.fa.fai' });
    expect(link).toHaveAttribute('href', `/v4/assets/${DIGEST}/files/hg38.fa.fai`);
  });

  it('saves a small file named by genome digest under the genome name', async () => {
    const GENOME = '5cb5b9d2d1e2b5e5e2f0a3c4d5e6f708';
    server.use(
      http.get(`${API}/assets/:digest`, () => HttpResponse.json(assetsFixture.results[1])),
      http.get(`${API}/assets/:digest/files`, () =>
        HttpResponse.json({ asset_digest: DIGEST, files: [`${GENOME}.chrom.sizes`] }),
      ),
      http.get(`${API}/assets/:digest/files/:file`, () => HttpResponse.text('chr1\t9\n')),
    );
    const saved: string[] = [];
    URL.createObjectURL = vi.fn(() => 'blob:sizes');
    URL.revokeObjectURL = vi.fn();
    const click = vi
      .spyOn(HTMLAnchorElement.prototype, 'click')
      .mockImplementation(function (this: HTMLAnchorElement) {
        saved.push(this.download);
      });
    renderAsset(serverConfig);
    await userEvent.click(
      await screen.findByRole('link', { name: 'Download hg38.chrom.sizes' }),
    );
    await waitFor(() => expect(saved).toEqual(['hg38.chrom.sizes']));
    click.mockRestore();
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
