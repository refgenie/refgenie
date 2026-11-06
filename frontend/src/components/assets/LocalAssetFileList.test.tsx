import { describe, expect, it } from 'vitest';
import { screen, waitFor, within } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { HttpResponse, http } from 'msw';
import { Route, Routes } from 'react-router-dom';
import { LocalAssetFileList } from './LocalAssetFileList';
import { AssetPage } from '../../pages/AssetPage';
import { server } from '../../test/server';
import { LOCAL_API } from '../../test/handlers';
import { localConfig, renderWithProviders, serverConfig } from '../../test/renderWithProviders';
import type { LocalAssetFiles, LocalFilePreview } from '../../types/api';
import type { UiConfig } from '../../types/ui';

const DIGEST = '1111111111111111111111111111111111111111111111111111111111111111';

const withLocalFiles: UiConfig = {
  ...localConfig,
  capabilities: { ...localConfig.capabilities, local_files: true },
};

const listing: LocalAssetFiles = {
  asset_digest: DIGEST,
  truncated: false,
  entries: [
    { path: 'hg38.fa.fai', type: 'file', size: 120, previewable: true },
    { path: 'hg38.fa.gz', type: 'file', size: 9000, previewable: false },
    { path: 'index', type: 'dir', size: null, previewable: false },
    { path: 'index/hg38.1.bt2', type: 'file', size: 4096, previewable: false },
  ],
};

function preview(overrides: Partial<LocalFilePreview> = {}): LocalFilePreview {
  return {
    path: 'hg38.fa.fai',
    size: 120,
    content: 'chr1\t248956422\t112\t70\t71\n',
    truncated: false,
    limit: 65536,
    ...overrides,
  };
}

function useListing(body: LocalAssetFiles = listing) {
  server.use(http.get(`${LOCAL_API}/files/assets/:digest`, () => HttpResponse.json(body)));
}

describe('LocalAssetFileList', () => {
  it('links every file to its local download route and lists directories', async () => {
    useListing();
    renderWithProviders(<LocalAssetFileList assetDigest={DIGEST} />, { config: withLocalFiles });

    const fai = await screen.findByRole('link', { name: 'Download hg38.fa.fai' });
    expect(fai).toHaveAttribute('href', `/v1/files/assets/${DIGEST}/raw/hg38.fa.fai`);
    // Nested paths keep their `/` separators.
    expect(screen.getByRole('link', { name: 'Download index/hg38.1.bt2' })).toHaveAttribute(
      'href',
      `/v1/files/assets/${DIGEST}/raw/index/hg38.1.bt2`,
    );
    // Directories are rows, not links.
    expect(screen.getByText('index/')).toBeInTheDocument();
    expect(screen.queryByRole('link', { name: /Download index$/ })).not.toBeInTheDocument();
  });

  it('offers View only for previewable files', async () => {
    useListing();
    renderWithProviders(<LocalAssetFileList assetDigest={DIGEST} />, { config: withLocalFiles });
    expect(await screen.findByRole('button', { name: 'View hg38.fa.fai' })).toBeInTheDocument();
    expect(screen.queryByRole('button', { name: 'View hg38.fa.gz' })).not.toBeInTheDocument();
    expect(screen.queryByRole('button', { name: 'View index/hg38.1.bt2' })).not.toBeInTheDocument();
  });

  it('opens a preview modal with the file contents', async () => {
    useListing();
    server.use(
      http.get(`${LOCAL_API}/files/assets/:digest/preview/*`, () => HttpResponse.json(preview())),
    );
    renderWithProviders(<LocalAssetFileList assetDigest={DIGEST} />, { config: withLocalFiles });

    await userEvent.click(await screen.findByRole('button', { name: 'View hg38.fa.fai' }));
    const contents = await screen.findByLabelText('Contents of hg38.fa.fai');
    expect(contents).toHaveTextContent('chr1');
    expect(screen.queryByText(/Truncated/)).not.toBeInTheDocument();

    await userEvent.click(screen.getAllByRole('button', { name: 'Close' })[0]);
    await waitFor(() =>
      expect(screen.queryByLabelText('Contents of hg38.fa.fai')).not.toBeInTheDocument(),
    );
  });

  it('says when a preview is truncated', async () => {
    useListing();
    server.use(
      http.get(`${LOCAL_API}/files/assets/:digest/preview/*`, () =>
        HttpResponse.json(preview({ truncated: true, size: 3 * 1024 * 1024 })),
      ),
    );
    renderWithProviders(<LocalAssetFileList assetDigest={DIGEST} />, { config: withLocalFiles });

    await userEvent.click(await screen.findByRole('button', { name: 'View hg38.fa.fai' }));
    expect(
      await screen.findByText(/Truncated: showing the first 64\.0 KB of 3\.0 MB/),
    ).toBeInTheDocument();
  });

  it('shows the server error when a preview is refused', async () => {
    useListing();
    server.use(
      http.get(`${LOCAL_API}/files/assets/:digest/preview/*`, () =>
        HttpResponse.json(
          {
            ok: false,
            error: { code: 'not_previewable', message: 'binary or compressed; download it' },
          },
          { status: 409 },
        ),
      ),
    );
    renderWithProviders(<LocalAssetFileList assetDigest={DIGEST} />, { config: withLocalFiles });

    await userEvent.click(await screen.findByRole('button', { name: 'View hg38.fa.fai' }));
    expect(await screen.findByRole('alert')).toHaveTextContent('download it');
  });
});

describe('AssetPage file section', () => {
  function renderAsset(config: UiConfig) {
    return renderWithProviders(
      <Routes>
        <Route path="/assets/:digest" element={<AssetPage />} />
      </Routes>,
      { route: `/assets/${DIGEST}`, config },
    );
  }

  it('reads files from disk when local_files is on', async () => {
    useListing();
    renderAsset(withLocalFiles);
    const section = (await screen.findByRole('heading', { name: 'Files' })).closest('section');
    expect(section).not.toBeNull();
    const link = await within(section as HTMLElement).findByRole('link', {
      name: 'Download hg38.fa.fai',
    });
    expect(link.getAttribute('href')).toMatch(/^\/v1\/files\//);
  });

  it('keeps the public server download links when local_files is off', async () => {
    renderAsset(serverConfig);
    const link = await screen.findByRole('link', { name: 'Download hg38.fa.fai' });
    expect(link).toHaveAttribute('href', `/v4/assets/${DIGEST}/files/hg38.fa.fai`);
  });
});
