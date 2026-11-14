import { describe, expect, it } from 'vitest';
import { render, screen, within } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { RouterProvider, createMemoryRouter } from 'react-router-dom';
import { LandingPage } from './LandingPage';
import { createRoutes } from '../app/routes';
import { ConfigProvider } from '../app/ConfigProvider';
import { ResourceCacheContext } from '../hooks/useResource';
import { createResourceCache } from '../services/resourceCache';
import { ToastProvider } from '../components/common/ToastProvider';
import {
  localConfig,
  renderWithProviders,
  serverConfig,
} from '../test/renderWithProviders';

/** The real router, so a hero search can actually land on `/genomes`. */
function renderApp() {
  const cache = createResourceCache();
  const router = createMemoryRouter(createRoutes(serverConfig), {
    initialEntries: ['/'],
  });
  return render(
    <ResourceCacheContext.Provider value={cache}>
      <ConfigProvider config={serverConfig}>
        <ToastProvider>
          <RouterProvider router={router} />
        </ToastProvider>
      </ConfigProvider>
    </ResourceCacheContext.Provider>,
  );
}

describe('LandingPage', () => {
  it('leads with what refgenie is', () => {
    renderWithProviders(<LandingPage />, { config: serverConfig });
    expect(screen.getByRole('heading', { level: 1 })).toHaveTextContent(
      /Reference genome data management/i,
    );
    expect(screen.getByText(/downloads, builds and organizes/i)).toBeInTheDocument();
  });

  it('shows the counts as one line once they load', async () => {
    renderWithProviders(<LandingPage />, { config: serverConfig });
    expect(await screen.findByText(/2 genomes · \d+ asset groups/)).toBeInTheDocument();
  });

  it('deep-links the hero search into the genome list', async () => {
    renderApp();
    await userEvent.type(screen.getByLabelText('Search genomes'), 'hg38{Enter}');
    expect(
      await screen.findByRole('heading', { level: 1, name: 'Genomes' }),
    ).toBeInTheDocument();
  });

  it('repeats none of the nav links', () => {
    renderWithProviders(<LandingPage />, { config: serverConfig });
    for (const name of [
      /Tree of life/,
      /^Genomes/,
      /^Species/,
      /^Recipes/,
      /^Aliases/,
    ]) {
      expect(screen.queryByRole('link', { name })).not.toBeInTheDocument();
    }
  });

  it('shows the same four install steps in both modes', () => {
    renderWithProviders(<LandingPage />, { config: localConfig });
    const steps = within(
      screen.getByRole('heading', { name: 'Get started' }).closest('section')!,
    ).getAllByRole('listitem');
    // The install line follows the version the server reports, and the
    // first pull is the small FASTA index, not the multi-GB FASTA.
    expect(steps.map((step) => step.querySelector('code')?.textContent)).toEqual([
      'pip install "refgenie>=1.0.0a2"',
      'refgenie init',
      'refgenie pull hg38/fasta_index',
      'refgenie seek hg38/fasta_index.chrom_sizes',
    ]);
    expect(screen.getByText(/About 12 KB/)).toBeInTheDocument();
  });

  it('says 1.0 is an alpha and links the details', () => {
    renderWithProviders(<LandingPage />, { config: serverConfig });
    expect(screen.getByText(/is an alpha release/)).toHaveTextContent(
      /pin the version in production pipelines/,
    );
    expect(screen.getByRole('link', { name: /Details/ })).toHaveAttribute(
      'href',
      'https://docs.refgenie.org/refgenie/',
    );
  });

  it('puts the AI section after Get started', () => {
    renderWithProviders(<LandingPage />, { config: serverConfig });
    const headings = screen.getAllByRole('heading', { level: 2 }).map((h) => h.textContent);
    expect(headings).toEqual(['Get started', 'Use it with an AI assistant']);
  });

  it('offers the AI section in both modes', () => {
    renderWithProviders(<LandingPage />, { config: serverConfig });
    expect(
      screen.getByRole('heading', { name: /Use it with an AI assistant/i }),
    ).toBeInTheDocument();
    expect(
      screen.getByText(/claude mcp add --transport http refgenie .*\/mcp$/),
    ).toBeInTheDocument();
  });

  it('offers the stdio MCP command on a local refgenie', () => {
    renderWithProviders(<LandingPage />, { config: localConfig });
    expect(
      screen.getByText('claude mcp add refgenie refgenie-mcp'),
    ).toBeInTheDocument();
  });

  it('links SKILL.md as a real file, not a client route', () => {
    renderWithProviders(<LandingPage />, { config: localConfig });
    const link = screen.getByRole('link', { name: 'SKILL.md' });
    expect(link).toHaveAttribute('href', '/SKILL.md');
    // No target=_blank: it is same-origin, and ExternalLink would be wrong here.
    expect(link).not.toHaveAttribute('target');
  });
});
