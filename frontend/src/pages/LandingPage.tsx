/**
 * The front door. Served at `/` in both modes.
 *
 * The stat tiles read `pagination.total` from the four core-router list
 * endpoints rather than `/v4/summary`, which the version4 router mounts only
 * in server mode (refgenie/server/main.py). One page, no mode branch: the
 * differences are capability-gated, never keyed off `config.mode`.
 */

import { useState } from 'react';
import type { FormEvent } from 'react';
import { Link, useNavigate } from 'react-router-dom';
import { useGenomes } from '../hooks/queries/useGenomes';
import { useAssetGroups } from '../hooks/queries/useAssetGroups';
import { useAssetClasses } from '../hooks/queries/useAssetClasses';
import { useRecipes } from '../hooks/queries/useRecipes';
import { useUiConfig } from '../hooks/useUiConfig';
import { useCapability } from '../hooks/useCapability';
import { useDocumentTitle } from '../hooks/useDocumentTitle';
import { CopyButton } from '../components/common/CopyButton';
import { ExternalLink } from '../components/common/ExternalLink';
import { Icon } from '../components/common/Icon';
import { ConnectCard } from '../components/bridge/ConnectCard';
import { PROJECT_LINKS } from '../services/projectLinks';
import type { CapabilityKey } from '../types/ui';

/** One count request each. `limit=1` is the smallest legal page. */
const COUNT_ONLY = { limit: 1 } as const;

const STEPS = [
  { label: 'Install', command: 'pip install refgenie', note: 'Python 3.11 or newer.' },
  {
    label: 'Pull an asset',
    command: 'refgenie pull hg38/fasta',
    note: 'Downloads and unpacks from the public server.',
  },
  {
    label: 'Use it',
    command: 'refgenie seek hg38/fasta',
    note: 'Prints the local path, ready to paste into a pipeline.',
  },
] as const;

/** Registers the read-only stdio MCP server with Claude Code. */
const MCP_COMMAND = 'claude mcp add refgenie refgenie-mcp';

interface Destination {
  to: string;
  title: string;
  body: string;
  /** Hidden when the capability is off. */
  capability?: CapabilityKey;
}

const DESTINATIONS: Destination[] = [
  {
    to: '/genomes',
    title: 'Genomes',
    body: 'Every assembly this server holds, with aliases, accessions, and the assets built for each one.',
  },
  {
    to: '/species',
    title: 'Species',
    body: 'Every organism represented here, by scientific and common name, with its genomes.',
  },
  {
    to: '/asset-classes',
    title: 'Asset classes',
    body: 'The typed shapes assets come in, the files each one contains, and the genomes that have them.',
  },
  {
    to: '/recipes',
    title: 'Recipes',
    body: 'The reproducible commands that build each asset class from its inputs.',
  },
  {
    to: '/tree',
    title: 'Tree of life',
    body: 'Browse the collection by taxonomy on a radial tree, then jump from a species to its genomes.',
    capability: 'archives',
  },
  {
    to: '/aliases',
    title: 'Aliases',
    body: 'The human-readable names that resolve to genome digests.',
  },
];

interface StatTileProps {
  value?: number;
  label: string;
  /** Omitted where the collection has no browse route of its own. */
  to?: string;
}

function StatTile({ value, label, to }: StatTileProps) {
  const body = (
    <>
      <span className="rg-stat__value">{value === undefined ? '—' : value.toLocaleString()}</span>
      <span className="rg-stat__label">{label}</span>
    </>
  );
  if (!to) return <div className="rg-stat">{body}</div>;
  return (
    <Link className="rg-stat rg-link--plain no-underline" to={to}>
      {body}
    </Link>
  );
}

export function LandingPage() {
  // The bare service name: `/` is the instance, not a subject inside it.
  useDocumentTitle(undefined);
  const config = useUiConfig();
  const navigate = useNavigate();
  const [query, setQuery] = useState('');
  const canSeeTree = useCapability('archives');

  const genomes = useGenomes(COUNT_ONLY);
  const assetGroups = useAssetGroups(COUNT_ONLY);
  const assetClasses = useAssetClasses(COUNT_ONLY);
  const recipes = useRecipes(COUNT_ONLY);

  const submit = (event: FormEvent) => {
    event.preventDefault();
    const q = query.trim();
    if (q !== '') navigate(`/genomes?q=${encodeURIComponent(q)}`);
  };

  const destinations = DESTINATIONS.filter(
    (d) => !d.capability || config.capabilities[d.capability],
  );

  // Root-relative like AboutPage's `${config.root_path}/docs`, so a sub-path
  // deployment resolves correctly. The absolute form is what a user pastes
  // into an assistant that can only fetch URLs.
  const skillPath = `${config.root_path}/SKILL.md`;
  const skillUrl = new URL(skillPath, window.location.origin).href;
  const skillPrompt = `Read ${skillUrl}, then tell me which genomes this refgenie server has.`;

  return (
    <div className="flex flex-col gap-12">
      {/* 1. Hero */}
      <section className="rg-hero text-center flex flex-col items-center gap-4">
        <img className="rg-hero__logo" src="refgenie_logo.svg" alt="" />
        <h1 className="text-4xl font-bold">Reference genome assets, ready to pull</h1>
        <p className="rg-hero__lede rg-muted text-lg">
          Refgenie is a standardized genome asset management system. It indexes, stores, and
          serves reference genome resources — aligner indexes, annotations, chromosome sizes —
          so every pipeline gets the same files, under the same names, wherever it runs.
        </p>

        <form className="rg-hero__search flex gap-2 w-full" onSubmit={submit}>
          <label className="sr-only" htmlFor="landing-search">
            Search genomes
          </label>
          <input
            id="landing-search"
            className="rg-field__input"
            type="search"
            placeholder="Search genomes — try hg38"
            value={query}
            onChange={(event) => setQuery(event.target.value)}
          />
          <button className="rg-btn rg-btn--primary" type="submit">
            Search
          </button>
        </form>

        <div className="flex flex-wrap justify-center gap-2">
          <Link className="rg-btn" to="/genomes">
            Browse genomes
          </Link>
          {canSeeTree && (
            <Link className="rg-btn" to="/tree">
              Tree of life
            </Link>
          )}
          <ExternalLink className="rg-btn" href={PROJECT_LINKS.docs}>
            Documentation
          </ExternalLink>
        </div>
      </section>

      {/* 2. Stat strip */}
      <section aria-label="What this instance holds">
        <div className="grid grid-cols-2 lg:grid-cols-4 gap-4">
          <StatTile value={genomes.data?.pagination.total} label="Genomes" to="/genomes" />
          <StatTile value={assetGroups.data?.pagination.total} label="Asset groups" />
          <StatTile
            value={assetClasses.data?.pagination.total}
            label="Asset classes"
            to="/asset-classes"
          />
          <StatTile value={recipes.data?.pagination.total} label="Recipes" to="/recipes" />
        </div>
      </section>

      {/* 3. Get started */}
      <section className="flex flex-col gap-4">
        <h2 className="text-xl font-semibold">Get started</h2>
        <ol className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {STEPS.map((step, index) => (
            <li className="rg-step flex flex-col gap-2" key={step.command}>
              <span className="flex items-center gap-2">
                <span className="rg-step__index">{index + 1}</span>
                <span className="font-semibold">{step.label}</span>
              </span>
              <span className="rg-command flex items-center gap-2">
                <code className="rg-code rg-code--inline flex-1 overflow-x-auto">
                  {step.command}
                </code>
                <CopyButton value={step.command} label={`Copy: ${step.command}`} />
              </span>
              <span className="rg-muted text-sm">{step.note}</span>
            </li>
          ))}
        </ol>
      </section>

      {/*
        4. Local refgenie. The card hides itself once connected, and on a local
        dash the bridge never mounts at all — this page IS the local refgenie.
      */}
      {config.mode !== 'local' && <ConnectCard />}

      {/*
        5. Use it with AI. Not gated on config.mode: the stdio MCP server and
        the SKILL.md URL are correct in both modes. The HTTP MCP endpoint is
        server-only, so it lives in SKILL.md, which spells out the mode check.
      */}
      <section className="flex flex-col gap-4">
        <h2 className="text-xl font-semibold">Use it with an AI assistant</h2>
        <p className="rg-muted text-sm">
          Refgenie is built to be driven by an AI assistant. Connect the read-only MCP server so
          an assistant can answer questions about your own genomes, or hand it the link below and
          it will learn this server&rsquo;s API on its own.
        </p>

        <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
          <div className="rg-card rg-card__body flex flex-col gap-2">
            <span className="font-semibold">Connect the MCP server</span>
            <span className="rg-muted text-sm">
              Read-only access to your local refgenie database. It can look things up; it never
              changes anything.
            </span>
            <span className="rg-command flex items-center gap-2">
              <code className="rg-code rg-code--inline flex-1 overflow-x-auto">{MCP_COMMAND}</code>
              <CopyButton value={MCP_COMMAND} label={`Copy: ${MCP_COMMAND}`} />
            </span>
          </div>

          <div className="rg-card rg-card__body flex flex-col gap-2">
            <span className="font-semibold">Point an assistant at this server</span>
            <span className="rg-muted text-sm">
              {/* A plain <a>, not <Link>: /SKILL.md is a static file on this
                  origin, so the SPA router must never see it, and it stays in
                  this tab. The `code` glyph marks it as a raw document this
                  instance serves rather than a page in this app; the visible
                  label already names the file, so it needs no sr-only text. */}
              <a className="rg-link" href={skillPath}>
                SKILL.md
                <Icon name="code" className="rg-icon--trailing" />
              </a>{' '}
              teaches any assistant this instance&rsquo;s API in one page. Paste this into a chat:
            </span>
            <span className="rg-command flex items-center gap-2">
              <code className="rg-code rg-code--inline flex-1 overflow-x-auto">{skillPrompt}</code>
              <CopyButton value={skillPrompt} label="Copy prompt" />
            </span>
          </div>
        </div>

        <p className="rg-muted text-sm">
          Prefer to drive the command-line tool?{' '}
          <ExternalLink href={PROJECT_LINKS.toolSkill}>The refgenie skill file</ExternalLink> covers
          install, pull, build and seek. See also{' '}
          <ExternalLink href={PROJECT_LINKS.useWithAi}>
            how to use refgenie with an AI assistant
          </ExternalLink>
          .
        </p>
      </section>

      {/* 6. Explore */}
      <section className="flex flex-col gap-4">
        <h2 className="text-xl font-semibold">Explore</h2>
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
          {destinations.map((d) => (
            <Link
              className="rg-card rg-card__body rg-link--plain no-underline flex flex-col gap-2"
              to={d.to}
              key={d.to}
            >
              <span className="font-semibold">{d.title}</span>
              <span className="rg-muted text-sm">{d.body}</span>
            </Link>
          ))}
        </div>
      </section>

      {/* 7. Learn more */}
      <section className="rg-muted text-sm">
        Running your own instance?{' '}
        <Link className="rg-link" to="/about">
          See what this one supports
        </Link>
        , or read the{' '}
        <ExternalLink href={PROJECT_LINKS.docs}>documentation</ExternalLink>.
      </section>
    </div>
  );
}
