import { describe, expect, it } from 'vitest';
import { screen, waitFor, within } from '@testing-library/react';
import { AboutPage } from './AboutPage';
import { localConfig, renderWithProviders, serverConfig } from '../test/renderWithProviders';

const project = () => screen.getByRole('region', { name: 'Refgenie, the project' });
const instance = () => screen.getByRole('region', { name: 'This instance' });

describe('AboutPage', () => {
  it('covers both concerns: the project and this instance', () => {
    renderWithProviders(<AboutPage />, { config: serverConfig });
    expect(screen.getByRole('heading', { level: 1 })).toHaveTextContent('About refgenie');
    expect(project()).toBeInTheDocument();
    expect(instance()).toBeInTheDocument();
  });

  it('sends readers to docs.refgenie.org, never the retired databio.org site', () => {
    renderWithProviders(<AboutPage />, { config: serverConfig });
    const hrefs = screen.getAllByRole('link').map((link) => link.getAttribute('href') ?? '');
    expect(hrefs.some((href) => href.startsWith('https://docs.refgenie.org'))).toBe(true);
    expect(hrefs.some((href) => href.includes('refgenie.databio.org'))).toBe(false);
  });

  it('never links the private refgenie1 repo', () => {
    renderWithProviders(<AboutPage />, { config: serverConfig });
    const hrefs = screen.getAllByRole('link').map((link) => link.getAttribute('href') ?? '');
    expect(hrefs.some((href) => href.includes('refgenie/refgenie1'))).toBe(false);
  });

  it('offers both citations with resolvable DOIs', () => {
    renderWithProviders(<AboutPage />, { config: serverConfig });
    expect(
      within(project()).getByRole('link', { name: /reference genome resource manager/i }),
    ).toHaveAttribute('href', 'https://doi.org/10.1093/gigascience/giz149');
    expect(
      within(project()).getByRole('link', { name: /Identity and compatibility/i }),
    ).toHaveAttribute('href', 'https://doi.org/10.1093/nargab/lqab036');
  });

  it('shows the project section unchanged on a local dash', () => {
    renderWithProviders(<AboutPage />, { config: localConfig });
    expect(
      within(project()).getByRole('heading', { name: 'How to cite refgenie' }),
    ).toBeInTheDocument();
    expect(
      within(project()).getByRole('link', { name: /Use the dashboard/ }),
    ).toBeInTheDocument();
  });

  it('names the instance by service name, not by a mode branch', () => {
    renderWithProviders(<AboutPage />, { config: localConfig });
    // The lede interpolates service_name. `refgenie local dashboard` also
    // appears in the Identity list, so match the sentence, not the name.
    expect(
      within(instance()).getByText(/Everything below describes the refgenie local dashboard/),
    ).toBeInTheDocument();
  });

  it('gates the standards endpoints on capabilities, not on mode', () => {
    renderWithProviders(<AboutPage />, { config: serverConfig });
    expect(
      within(instance()).getByRole('link', { name: /Sequence collections service info/ }),
    ).toBeInTheDocument();
    expect(
      within(instance()).getByRole('link', { name: /GA4GH DRS service info/ }),
    ).toBeInTheDocument();
  });

  it('hides them on a dash, where neither router is mounted', () => {
    renderWithProviders(<AboutPage />, { config: localConfig });
    expect(
      within(instance()).queryByRole('link', { name: /Sequence collections service info/ }),
    ).not.toBeInTheDocument();
    expect(
      within(instance()).queryByRole('link', { name: /GA4GH DRS service info/ }),
    ).not.toBeInTheDocument();
  });

  it('points endpoints at the API origin on a cross-origin deployment', async () => {
    renderWithProviders(<AboutPage />, {
      config: { ...serverConfig, api_base: 'https://api.refgenie.org/v4' },
    });
    await waitFor(() =>
      expect(
        within(instance()).getByRole('link', { name: /OpenAPI schema/ }),
      ).toHaveAttribute('href', 'https://api.refgenie.org/openapi.json'),
    );
  });

  it('still lists every capability flag', () => {
    renderWithProviders(<AboutPage />, { config: serverConfig });
    expect(within(instance()).getByText('aliases_write')).toBeInTheDocument();
    expect(within(instance()).getByText('drs')).toBeInTheDocument();
  });

  it('links SKILL.md as a real file, not a client route', () => {
    renderWithProviders(<AboutPage />, { config: localConfig });
    const link = within(instance()).getByRole('link', { name: /SKILL\.md/ });
    expect(link).toHaveAttribute('href', '/SKILL.md');
    expect(link).not.toHaveAttribute('target');
  });
});
