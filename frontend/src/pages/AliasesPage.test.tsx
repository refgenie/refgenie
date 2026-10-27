/**
 * The Aliases page lists aliases everywhere and adds them where the instance
 * allows: the button and its modal are gated on `aliases_write`, so the public
 * site keeps a read-only list.
 */

import { screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { HttpResponse, http } from 'msw';
import { describe, expect, it } from 'vitest';
import { AliasesPage } from './AliasesPage';
import { renderWithProviders, serverConfig } from '../test/renderWithProviders';
import { server } from '../test/server';
import { API, LOCAL_API } from '../test/handlers';

function seedGenomes() {
  server.use(
    http.get(`${API}/genomes`, () =>
      HttpResponse.json({
        items: [
          {
            digest: 'genome-digest-1',
            aliases: ['hg38'],
            description: null,
            asset_count: 1,
            species_name: null,
            common_name: null,
            taxon_id: null,
            assembly_source: null,
            assembly_accession: null,
          },
        ],
        pagination: { offset: 0, limit: 25, total: 1 },
      }),
    ),
  );
}

describe('AliasesPage', () => {
  it('adds an alias from the page head with the action header', async () => {
    seedGenomes();
    let headerSeen: string | null = null;
    let body: unknown;
    server.use(
      http.post(`${LOCAL_API}/actions/aliases`, async ({ request }) => {
        headerSeen = request.headers.get('x-refgenie-action');
        body = await request.json();
        return HttpResponse.json({ ok: true, message: 'Alias set.', data: null });
      }),
    );

    renderWithProviders(<AliasesPage />, { route: '/aliases' });
    await userEvent.click(await screen.findByRole('button', { name: 'Add alias' }));

    const heading = await screen.findByRole('heading', { name: 'Add alias' });
    await userEvent.type(screen.getByLabelText(/^Alias \*/), 'GRCh38');
    await userEvent.selectOptions(await screen.findByLabelText(/^Genome/), 'hg38');
    // Two "Add alias" buttons now: the page head's and the modal's submit.
    const submit = screen.getAllByRole('button', { name: 'Add alias' }).find(
      (button) => button.getAttribute('type') === 'submit',
    );
    await userEvent.click(submit!);

    await waitFor(() => expect(body).toEqual({ alias: 'GRCh38', genome_digest: 'genome-digest-1' }));
    expect(headerSeen).not.toBeNull();
    await waitFor(() => expect(heading).not.toBeInTheDocument());
  });

  it('offers no Add alias on an instance without aliases_write', async () => {
    renderWithProviders(<AliasesPage />, { route: '/aliases', config: serverConfig });
    await screen.findByRole('table');
    expect(screen.queryByRole('button', { name: 'Add alias' })).not.toBeInTheDocument();
  });
});
