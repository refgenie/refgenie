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
import { aliasesFixture } from '../test/fixtures';

function seedGenomes() {
  server.use(
    http.get(`${API}/genomes`, () =>
      HttpResponse.json({
        results: [
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
            species: null,
          },
        ],
        pagination: { page: 0, page_size: 25, total: 1 },
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
    const submit = screen
      .getAllByRole('button', { name: 'Add alias' })
      .find((button) => button.getAttribute('type') === 'submit');
    await userEvent.click(submit!);

    await waitFor(() =>
      expect(body).toEqual({ alias: 'GRCh38', genome_digest: 'genome-digest-1' }),
    );
    expect(headerSeen).not.toBeNull();
    await waitFor(() => expect(heading).not.toBeInTheDocument());
  });

  it('reads a 1-based page from the URL and falls back to the first page past the end', async () => {
    const requested: Array<string | null> = [];
    server.use(
      http.get(`${API}/aliases`, ({ request }) => {
        const page = new URL(request.url).searchParams.get('page');
        requested.push(page);
        // The server's answer for `page * page_size >= total`.
        if (Number(page) > 0) {
          return HttpResponse.json({ detail: 'page is past the end' }, { status: 400 });
        }
        return HttpResponse.json(aliasesFixture);
      }),
    );

    renderWithProviders(<AliasesPage />, {
      route: '/aliases?page=3',
      config: serverConfig,
    });
    expect(await screen.findByRole('link', { name: 'hg38' })).toBeInTheDocument();
    // `?page=3` is API page 2; the 400 sends the view back to API page 0.
    expect(requested).toEqual(['2', '0']);
  });

  it('offers no Add alias on an instance without aliases_write', async () => {
    renderWithProviders(<AliasesPage />, { route: '/aliases', config: serverConfig });
    await screen.findByRole('table');
    expect(screen.queryByRole('button', { name: 'Add alias' })).not.toBeInTheDocument();
  });

  it("shows each alias's species, assembly and asset count from the genome list", async () => {
    renderWithProviders(<AliasesPage />, { route: '/aliases', config: serverConfig });
    const row = (await screen.findByRole('link', { name: 'hg38' })).closest('tr')!;
    await waitFor(() => expect(row).toHaveTextContent('Homo sapiens'));
    expect(row).toHaveTextContent('NCBI · GCA_000001405.15');
    expect(row).toHaveTextContent('3');
  });
});
