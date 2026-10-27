import { describe, expect, it } from 'vitest';
import { render, screen, waitFor, within } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { RouterProvider, createMemoryRouter } from 'react-router-dom';
import { LandingPage } from './LandingPage';
import { createRoutes } from '../app/routes';
import { ConfigProvider } from '../app/ConfigProvider';
import { ResourceCacheContext } from '../services/resourceCacheContext';
import { createResourceCache, NO_RETRY } from '../services/resourceCache';
import { ToastProvider } from '../components/common/ToastProvider';
import { localConfig, renderWithProviders, serverConfig } from '../test/renderWithProviders';

/** The stat strip, which is the only region on the page. */
const stats = () => screen.getByRole('region', { name: 'What this instance holds' });

/** The real router, so a hero search can actually land on `/genomes`. */
function renderApp() {
  const cache = createResourceCache({ retry: NO_RETRY });
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
    expect(screen.getByRole('heading', { level: 1 })).toHaveTextContent(/ready to pull/i);
    expect(
      screen.getByText(/standardized genome asset management system/i),
    ).toBeInTheDocument();
  });

  it('fills the stat tiles from pagination totals', async () => {
    renderWithProviders(<LandingPage />, { config: serverConfig });
    const genomes = within(stats()).getByRole('link', { name: /Genomes/ });
    await waitFor(() => expect(genomes).not.toHaveTextContent('—'));
    expect(genomes).toHaveTextContent('2');
  });

  it('deep-links the hero search into the genome list', async () => {
    renderApp();
    await userEvent.type(screen.getByLabelText('Search genomes'), 'hg38{Enter}');
    expect(
      await screen.findByRole('heading', { level: 1, name: 'Genomes' }),
    ).toBeInTheDocument();
  });

  it('hides the tree in local mode, where the species summary does not exist', () => {
    renderWithProviders(<LandingPage />, { config: localConfig });
    expect(screen.queryByRole('link', { name: /Tree of life/ })).not.toBeInTheDocument();
  });

  it('offers the tree in server mode', () => {
    renderWithProviders(<LandingPage />, { config: serverConfig });
    expect(screen.getAllByRole('link', { name: /Tree of life/ }).length).toBeGreaterThan(0);
  });

  it('shows the same install steps in both modes', () => {
    renderWithProviders(<LandingPage />, { config: localConfig });
    expect(screen.getByText('pip install refgenie')).toBeInTheDocument();
    expect(screen.getByText('refgenie pull hg38/fasta')).toBeInTheDocument();
  });

  it('offers the AI section in both modes', () => {
    renderWithProviders(<LandingPage />, { config: serverConfig });
    expect(
      screen.getByRole('heading', { name: /Use it with an AI assistant/i }),
    ).toBeInTheDocument();
    expect(screen.getByText('claude mcp add refgenie refgenie-mcp')).toBeInTheDocument();
  });

  it('links SKILL.md as a real file, not a client route', () => {
    renderWithProviders(<LandingPage />, { config: localConfig });
    const link = screen.getByRole('link', { name: 'SKILL.md' });
    expect(link).toHaveAttribute('href', '/SKILL.md');
    // No target=_blank: it is same-origin, and ExternalLink would be wrong here.
    expect(link).not.toHaveAttribute('target');
  });

  it('has no flat Assets destination, because /assets is not a browse route', () => {
    renderWithProviders(<LandingPage />, { config: serverConfig });
    expect(screen.queryByRole('link', { name: /^Assets/ })).not.toBeInTheDocument();
  });
});
