import { describe, expect, it } from 'vitest';
import { screen } from '@testing-library/react';
import { HttpResponse, http } from 'msw';
import { CollectionsPage } from './CollectionsPage';
import { server } from '../test/server';
import { API } from '../test/handlers';
import { renderWithProviders, serverConfig } from '../test/renderWithProviders';
import { vgpCollection } from '../test/collectionFixtures';

describe('CollectionsPage', () => {
  it('shows one card per collection with its genome count', async () => {
    server.use(
      http.get(`${API}/genomes`, ({ request }) => {
        const params = new URL(request.url).searchParams;
        const total = params.get('store') === 'vgp' ? 605 : 0;
        return HttpResponse.json({
          results: [],
          pagination: { page: 0, page_size: 1, total },
        });
      }),
    );
    const other = {
      ...vgpCollection,
      id: 'hprc',
      name: 'HPRC assemblies',
      stores: ['hprc'],
    };
    renderWithProviders(<CollectionsPage />, {
      config: { ...serverConfig, collections: [vgpCollection, other] },
    });
    expect(screen.getAllByRole('listitem')).toHaveLength(2);
    expect(screen.getByRole('link', { name: vgpCollection.name })).toHaveAttribute(
      'href',
      '/collections/vgp',
    );
    expect(await screen.findByText('605 genomes')).toBeInTheDocument();
    expect(await screen.findByText('0 genomes')).toBeInTheDocument();
  });

  it('shows an empty state with no collections', () => {
    renderWithProviders(<CollectionsPage />, { config: serverConfig });
    expect(
      screen.getByText('This server has no collections configured.'),
    ).toBeInTheDocument();
  });
});
