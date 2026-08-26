/**
 * The single nav registry. Adding a managed surface later is one line here
 * plus a route; NavBar filters on capability, never on mode.
 */

import type { CapabilityKey } from '../types/ui';

export interface NavEntry {
  to: string;
  label: string;
  capability?: CapabilityKey;
}

export const NAV: NavEntry[] = [
  { to: '/genomes', label: 'Genomes' },
  { to: '/assets', label: 'Assets' },
  { to: '/asset-classes', label: 'Asset classes' },
  { to: '/recipes', label: 'Recipes' },
  { to: '/aliases', label: 'Aliases' },
  { to: '/remote', label: 'Remote', capability: 'remote_browse' },
  { to: '/build', label: 'Build', capability: 'build' },
  { to: '/manage', label: 'Manage', capability: 'subscriptions' },
  { to: '/jobs', label: 'Jobs', capability: 'jobs' },
  { to: '/about', label: 'About' },
];
