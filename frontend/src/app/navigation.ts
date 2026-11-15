/** The nav registry. NavBar filters entries by capability. */

import type { CapabilityKey, UiConfig } from '../types/ui';

export interface NavEntry {
  to: string;
  label: string;
  /** Shown on hover: what the page lists, in plain words. */
  hint?: string;
  capability?: CapabilityKey;
  /** Also listed, under `/local`, when browsing a connected local refgenie. */
  local?: boolean;
  /** Shown only when this returns true for the bootstrap config. */
  visible?: (config: UiConfig) => boolean;
}

export const NAV: NavEntry[] = [
  { to: '/', label: 'Home' },
  { to: '/genomes', label: 'Genomes', local: true },
  { to: '/species', label: 'Species', local: true },
  {
    to: '/collections',
    label: 'Collections',
    hint: 'Hand-picked sets of genomes, such as the VGP assemblies',
    visible: (c) => c.collections.length > 0,
  },
  {
    to: '/asset-classes',
    label: 'Asset classes',
    hint: 'The kinds of files refgenie manages, such as fasta or bowtie2_index',
    local: true,
  },
  {
    to: '/recipes',
    label: 'Recipes',
    hint: 'The build steps that make each kind of file',
    local: true,
  },
  {
    to: '/aliases',
    label: 'Aliases',
    hint: 'The names genomes go by, such as hg38',
    local: true,
  },
  { to: '/remote', label: 'Remote', capability: 'remote_browse' },
  { to: '/servers', label: 'Servers', capability: 'remote_browse', local: true },
  { to: '/build', label: 'Build', capability: 'build' },
  { to: '/manage', label: 'Manage', capability: 'subscriptions' },
  { to: '/jobs', label: 'Jobs', capability: 'jobs', local: true },
];
