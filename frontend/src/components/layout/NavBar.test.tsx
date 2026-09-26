import { describe, expect, it } from 'vitest';
import { screen } from '@testing-library/react';
import { NavBar } from './NavBar';
import { localConfig, renderWithProviders, serverConfig } from '../../test/renderWithProviders';

describe('NavBar', () => {
  it('always shows the browse entries', () => {
    renderWithProviders(<NavBar />, { config: serverConfig });
    for (const label of [
      'Home',
      'Genomes',
      'Species',
      'Asset classes',
      'Recipes',
      'Aliases',
      'About',
    ]) {
      expect(screen.getByRole('link', { name: label })).toBeInTheDocument();
    }
  });

  it('shows Species on a local dash too — it is not capability-gated', () => {
    // The species table is grouped from /v4/genomes, which both modes serve,
    // so unlike the tree it needs no `archives` capability.
    renderWithProviders(<NavBar />, { config: localConfig });
    expect(screen.getByRole('link', { name: 'Species' })).toBeInTheDocument();
  });

  it('has no flat Assets browse entry', () => {
    // `/assets` was deleted: asset classes are the browse entry point, and
    // `/assets/:digest` survives only as the detail route.
    renderWithProviders(<NavBar />, { config: serverConfig });
    expect(screen.queryByRole('link', { name: 'Assets' })).not.toBeInTheDocument();
  });

  it('shows the tree in server mode, where the archives capability is on', () => {
    renderWithProviders(<NavBar />, { config: serverConfig });
    expect(screen.getByRole('link', { name: 'Tree of life' })).toBeInTheDocument();
  });

  it('hides the tree in local mode, where /species answers the same question', () => {
    renderWithProviders(<NavBar />, { config: localConfig });
    expect(screen.queryByRole('link', { name: 'Tree of life' })).not.toBeInTheDocument();
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
