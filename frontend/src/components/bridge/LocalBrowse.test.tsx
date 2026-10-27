/**
 * The read-only `/local` views beyond genomes and assets: the same pages,
 * rendered inside `LocalScope`, reading the connected local refgenie and
 * keeping their links inside the scope.
 */
import { beforeEach, describe, expect, it } from 'vitest';
import { screen, within } from '@testing-library/react';
import { HttpResponse, http } from 'msw';
import { LocalScope } from './LocalScope';
import { ConnectCard } from './ConnectCard';
import { AliasesPage } from '../../pages/AliasesPage';
import { RecipesPage } from '../../pages/RecipesPage';
import { JobsPage } from '../../pages/JobsPage';
import { server } from '../../test/server';
import { aliasesFixture, recipesFixture } from '../../test/fixtures';
import { makeJob, page } from '../../test/jobFixtures';
import { renderWithProviders, serverConfig } from '../../test/renderWithProviders';
import { useBridgeStore } from '../../stores/bridgeStore';
import { coerceCapabilities } from '../../services/capabilities';
import { CAPABILITY_KEYS } from '../../types/ui';
import type { PingResponse } from '../../services/bridge/contract';

const BRIDGE_ORIGIN = 'http://localhost:8080';

/** A local instance that claims every capability, write ones included. */
const ping: PingResponse = {
  service: 'refgenie',
  bridge_version: 1,
  mode: 'local',
  refgenie_version: '1.0.0',
  instance_id: 'instance-1',
  instance_label: 'my laptop refgenie',
  bridge_mode: 'full',
  action_header: 'X-Refgenie-Action',
  capabilities: coerceCapabilities(Object.fromEntries(CAPABILITY_KEYS.map((key) => [key, true]))),
  bridge: { actions_cross_origin: true },
  counts: { genomes: 3, assets: 12 },
};

beforeEach(() => {
  localStorage.clear();
  useBridgeStore.getState().reset();
  useBridgeStore.getState().setConnected(ping, 8080, BRIDGE_ORIGIN);
});

describe('/local browse views', () => {
  it('lists local aliases with scoped links and no write controls', async () => {
    server.use(
      http.get(`${BRIDGE_ORIGIN}/v4/aliases`, () => HttpResponse.json(aliasesFixture)),
    );
    renderWithProviders(
      <LocalScope>
        <AliasesPage />
      </LocalScope>,
      { config: serverConfig, route: '/local/aliases' },
    );

    const link = await screen.findByRole('link', { name: 'hg38' });
    expect(link).toHaveAttribute('href', `/local/genomes/${aliasesFixture.items[0].genome_digest}`);
    expect(screen.queryByRole('button', { name: /add alias/i })).toBeNull();
  });

  it('lists local recipes with links that stay under /local', async () => {
    server.use(
      http.get(`${BRIDGE_ORIGIN}/v4/recipes`, () => HttpResponse.json(recipesFixture)),
    );
    renderWithProviders(
      <LocalScope>
        <RecipesPage />
      </LocalScope>,
      { config: serverConfig, route: '/local/recipes' },
    );
    const link = await screen.findByRole('link', { name: /fasta_txome/ });
    expect(link).toHaveAttribute('href', '/local/recipes/17');
  });

  it('shows the local job history read-only', async () => {
    server.use(
      http.get(`${BRIDGE_ORIGIN}/v1/jobs`, () =>
        HttpResponse.json(
          page([
            makeJob({
              status: 'succeeded',
              result: { asset_digest: 'asset-digest-1', registry_path: 'hg38/fasta:default' },
            }),
          ]),
        ),
      ),
    );
    renderWithProviders(
      <LocalScope>
        <JobsPage />
      </LocalScope>,
      { config: serverConfig, route: '/local/jobs' },
    );
    const result = await screen.findByRole('link', { name: 'hg38/fasta:default' });
    expect(result).toHaveAttribute('href', '/local/assets/asset-digest-1');
    expect(screen.queryByRole('button', { name: /cancel/i })).toBeNull();
  });

  it('offers local navigation between the scoped views', () => {
    renderWithProviders(
      <LocalScope>
        <p>content</p>
      </LocalScope>,
      { config: serverConfig, route: '/local/recipes' },
    );
    const nav = screen.getByRole('navigation', { name: /local refgenie/i });
    const hrefs = within(nav)
      .getAllByRole('link')
      .map((link) => link.getAttribute('href'));
    expect(hrefs).toEqual([
      '/local/genomes',
      '/local/species',
      '/local/asset-classes',
      '/local/recipes',
      '/local/aliases',
      '/local/jobs',
    ]);
    expect(within(nav).getByRole('link', { name: 'Recipes' })).toHaveClass(
      'rg-nav__item--active',
    );
  });

  it('turns the landing card into a way into /local once connected', () => {
    renderWithProviders(<ConnectCard />, { config: serverConfig });
    expect(screen.getByText(/3 genomes and 12 assets/)).toBeInTheDocument();
    expect(screen.getByRole('link', { name: /browse local genomes/i })).toHaveAttribute(
      'href',
      '/local/genomes',
    );
    expect(screen.getByRole('link', { name: /local jobs/i })).toHaveAttribute(
      'href',
      '/local/jobs',
    );
  });
});
