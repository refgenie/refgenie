/**
 * jsdom has no layout: `getBoundingClientRect()` returns zeros, so `drawTree`
 * bails and the SVG stays empty. That is correct behaviour, not a failure —
 * every bug this port fixes lives in the chrome, so the chrome is what is
 * asserted here.
 */

import { beforeAll, describe, expect, it } from 'vitest';
import { screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { TreeExplorerPage } from './TreeExplorerPage';
import { renderWithProviders, serverConfig } from '../test/renderWithProviders';

beforeAll(() => {
  // jsdom does not implement ResizeObserver, and the tree observes its panel.
  globalThis.ResizeObserver = class {
    observe() {}
    unobserve() {}
    disconnect() {}
  } as unknown as typeof ResizeObserver;
});

describe('TreeExplorerPage', () => {
  it('shows every control without a hover', () => {
    renderWithProviders(<TreeExplorerPage />, { config: serverConfig });
    for (const name of ['Zoom in', 'Zoom out', 'Reset', 'Inspect', 'Pan', 'Full screen']) {
      expect(screen.getByRole('button', { name })).toBeVisible();
    }
    expect(screen.getByLabelText('Filter species')).toBeVisible();
    expect(screen.getByLabelText('Group by taxonomic level')).toBeVisible();
  });

  it('starts in inspect mode, and the toggle is labelled', async () => {
    renderWithProviders(<TreeExplorerPage />, { config: serverConfig });
    expect(screen.getByRole('button', { name: 'Inspect' })).toHaveAttribute(
      'aria-pressed',
      'true',
    );
    await userEvent.click(screen.getByRole('button', { name: 'Pan' }));
    expect(screen.getByRole('button', { name: 'Pan' })).toHaveAttribute('aria-pressed', 'true');
    expect(screen.getByRole('button', { name: 'Inspect' })).toHaveAttribute(
      'aria-pressed',
      'false',
    );
  });

  it('does not enter full screen when you type', async () => {
    renderWithProviders(<TreeExplorerPage />, { config: serverConfig });
    await userEvent.type(screen.getByLabelText('Filter species'), 'Homo');
    expect(screen.getByRole('button', { name: 'Full screen' })).toBeInTheDocument();
    expect(screen.getByRole('heading', { level: 1, name: 'Tree of life' })).toBeInTheDocument();
  });

  it('leaves full screen on Escape', async () => {
    renderWithProviders(<TreeExplorerPage />, { config: serverConfig });
    await userEvent.click(screen.getByRole('button', { name: 'Full screen' }));
    expect(screen.getByRole('button', { name: 'Exit full screen' })).toBeInTheDocument();
    await userEvent.keyboard('{Escape}');
    expect(screen.getByRole('button', { name: 'Full screen' })).toBeInTheDocument();
  });

  it('selects a species from the URL despite the taxonomy spelling it differently', async () => {
    // taxa.json says "Homo Sapiens"; the API says "Homo sapiens". Before the
    // species index this lookup missed and the panel claimed there were no
    // human genomes.
    renderWithProviders(<TreeExplorerPage />, {
      config: serverConfig,
      route: '/tree?species=Homo%20sapiens',
    });
    expect(
      await screen.findByRole('heading', { level: 2, name: 'Homo Sapiens' }),
    ).toBeInTheDocument();
    expect(await screen.findByText(/2 genomes · 4 assets/)).toBeInTheDocument();
    expect(screen.getByRole('link', { name: 'Species table' })).toHaveAttribute(
      'href',
      '/species?q=Homo%20sapiens',
    );
  });

  it('links out to the genome list filtered by the species column', async () => {
    renderWithProviders(<TreeExplorerPage />, {
      config: serverConfig,
      route: '/tree?species=Homo%20sapiens',
    });
    expect(await screen.findByRole('link', { name: 'View genomes' })).toHaveAttribute(
      'href',
      '/genomes?q=Homo+sapiens&fields=species_name&op=contains',
    );
  });

  it('selects nothing for a species the taxonomy does not carry', async () => {
    renderWithProviders(<TreeExplorerPage />, {
      config: serverConfig,
      route: '/tree?species=Nonexistent%20species',
    });
    await screen.findByRole('heading', { level: 1, name: 'Tree of life' });
    expect(screen.queryByRole('heading', { level: 2 })).not.toBeInTheDocument();
  });

  it('offers all eight taxonomic levels, grouping by class', () => {
    renderWithProviders(<TreeExplorerPage />, { config: serverConfig });
    const select = screen.getByLabelText<HTMLSelectElement>('Group by taxonomic level');
    expect(select.value).toBe('class');
    expect(select.options).toHaveLength(8);
  });
});
