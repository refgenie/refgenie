import { describe, expect, it } from 'vitest';
import { screen } from '@testing-library/react';
import { NavBar } from './NavBar';
import { localConfig, renderWithProviders, serverConfig } from '../../test/renderWithProviders';

describe('NavBar', () => {
  it('always shows the browse entries', () => {
    renderWithProviders(<NavBar />, { config: serverConfig });
    for (const label of ['Genomes', 'Assets', 'Asset classes', 'Recipes', 'Aliases', 'About']) {
      expect(screen.getByRole('link', { name: label })).toBeInTheDocument();
    }
  });

  it('hides Remote, Manage and Jobs in server mode', () => {
    renderWithProviders(<NavBar />, { config: serverConfig });
    expect(screen.queryByRole('link', { name: 'Remote' })).not.toBeInTheDocument();
    expect(screen.queryByRole('link', { name: 'Manage' })).not.toBeInTheDocument();
    expect(screen.queryByRole('link', { name: 'Jobs' })).not.toBeInTheDocument();
  });

  it('shows Remote in local mode, where remote_browse is on', () => {
    renderWithProviders(<NavBar />, { config: localConfig });
    expect(screen.getByRole('link', { name: 'Remote' })).toBeInTheDocument();
  });
});
