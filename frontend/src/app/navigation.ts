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
  { to: '/', label: 'Home' },
  { to: '/genomes', label: 'Genomes' },
  // No capability: the species table is grouped from `/v4/genomes`, which both
  // modes serve, so a local dash gets it too.
  { to: '/species', label: 'Species' },
  { to: '/asset-classes', label: 'Asset classes' },
  { to: '/recipes', label: 'Recipes' },
  { to: '/aliases', label: 'Aliases' },
  // `archives` is `not is_local` (refgenie/server/main.py). A 560-taxon radial
  // tree over one user's handful of local genomes is noise; the flat /species
  // table is the local-mode answer to the same question. The route itself stays
  // registered in both modes.
  { to: '/tree', label: 'Tree of life', capability: 'archives' },
  { to: '/remote', label: 'Remote', capability: 'remote_browse' },
  { to: '/build', label: 'Build', capability: 'build' },
  { to: '/manage', label: 'Manage', capability: 'subscriptions' },
  { to: '/jobs', label: 'Jobs', capability: 'jobs' },
  { to: '/about', label: 'About' },
];
