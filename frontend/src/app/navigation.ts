/** The nav registry. NavBar filters entries by capability. */

import type { CapabilityKey } from '../types/ui';

export interface NavEntry {
  to: string;
  label: string;
  capability?: CapabilityKey;
  /** Also listed, under `/local`, when browsing a connected local refgenie. */
  local?: boolean;
}

export const NAV: NavEntry[] = [
  { to: '/', label: 'Home' },
  { to: '/genomes', label: 'Genomes', local: true },
  { to: '/species', label: 'Species', local: true },
  { to: '/asset-classes', label: 'Asset classes', local: true },
  { to: '/recipes', label: 'Recipes', local: true },
  { to: '/aliases', label: 'Aliases', local: true },
  { to: '/remote', label: 'Remote', capability: 'remote_browse' },
  { to: '/servers', label: 'Servers', capability: 'remote_browse', local: true },
  { to: '/build', label: 'Build', capability: 'build' },
  { to: '/manage', label: 'Manage', capability: 'subscriptions' },
  { to: '/jobs', label: 'Jobs', capability: 'jobs', local: true },
];
