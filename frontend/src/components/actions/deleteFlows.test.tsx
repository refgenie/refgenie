/**
 * The destructive paths. Both go through `ConfirmModal`; the genome one is
 * gated on typing the alias, because `genome.remove` cascades unconditionally.
 */

import { screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { HttpResponse, http } from 'msw';
import { describe, expect, it } from 'vitest';
import { DeleteAssetButton } from './DeleteAssetButton';
import { DeleteGenomeButton } from './DeleteGenomeButton';
import { renderWithProviders } from '../../test/renderWithProviders';
import { server } from '../../test/server';
import { LOCAL_API } from '../../test/handlers';

describe('DeleteAssetButton', () => {
  it('issues a DELETE carrying the action header and names the seek keys', async () => {
    let method: string | null = null;
    let headerSeen: string | null = null;
    server.use(
      http.delete(`${LOCAL_API}/actions/assets/:digest`, ({ request }) => {
        method = request.method;
        headerSeen = request.headers.get('x-refgenie-action');
        return HttpResponse.json({ digest: 'asset-1' });
      }),
    );

    renderWithProviders(
      <DeleteAssetButton
        digest="asset-1"
        registryPath="fasta:default"
        seekKeys={[
          { name: 'fasta', value: '/x.fa', description: null, type: 'file', size: null },
        ]}
      />,
    );

    await userEvent.click(screen.getByRole('button', { name: 'Delete' }));
    expect(screen.getByText('fasta:default')).toBeInTheDocument();
    expect(screen.getByText(/seek keys will disappear/)).toBeInTheDocument();

    await userEvent.click(screen.getByRole('button', { name: 'Delete asset' }));
    await waitFor(() => expect(method).toBe('DELETE'));
    expect(headerSeen).not.toBeNull();
  });
});

describe('DeleteGenomeButton', () => {
  it('keeps confirm disabled until the alias is typed, then deletes', async () => {
    let called = false;
    server.use(
      http.delete(`${LOCAL_API}/actions/genomes/:ref`, () => {
        called = true;
        return HttpResponse.json({ digest: 'genome-1', removed_assets: 3 });
      }),
    );

    renderWithProviders(
      <DeleteGenomeButton digest="genome-1" primaryAlias="hg38" aliasCount={2} assetCount={3} />,
    );

    await userEvent.click(screen.getByRole('button', { name: 'Delete genome' }));
    // The modal states the real consequence, with counts.
    expect(screen.getByText(/all 3 assets/)).toBeInTheDocument();

    const confirm = screen.getByRole('button', {
      name: 'Delete genome and all its assets',
    });
    expect(confirm).toBeDisabled();

    await userEvent.type(screen.getByLabelText(/to confirm/), 'hg38');
    expect(confirm).toBeEnabled();

    await userEvent.click(confirm);
    await waitFor(() => expect(called).toBe(true));
  });
});
