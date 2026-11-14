import { describe, expect, it } from 'vitest';
import { render, screen } from '@testing-library/react';
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

  it('shows the same install steps in both modes', () => {
    renderWithProviders(<LandingPage />, { config: localConfig });
    // The install line follows the version the server reports.
    expect(screen.getByText('pip install "refgenie>=1.0.0a2"')).toBeInTheDocument();
    expect(screen.getByText('refgenie pull hg38/fasta')).toBeInTheDocument();
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
