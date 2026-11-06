/**
 * The landing page at `/`. The counts read `pagination.total` from list
 * endpoints, not `/v4/summary`, which only server mode has.
 */

import { useState } from 'react';
import type { FormEvent } from 'react';
import { useNavigate } from 'react-router-dom';
import { useGenomes } from '../hooks/queries/useGenomes';
import { useAssetGroups } from '../hooks/queries/useAssetGroups';
import { useAssetClasses } from '../hooks/queries/useAssetClasses';
import { useRecipes } from '../hooks/queries/useRecipes';
import { useUiConfig } from '../hooks/useUiConfig';
import { useApiClient } from '../hooks/useApiClient';
import { useDocumentTitle } from '../hooks/useDocumentTitle';
import { CopyButton } from '../components/common/CopyButton';
import { ExternalLink } from '../components/common/ExternalLink';
import { Icon } from '../components/common/Icon';
import { PROJECT_LINKS } from '../services/projectLinks';
import { absoluteServerRoot } from '../services/instanceLinks';

/** Count-only request: one row per page. */
const COUNT_ONLY = { pageSize: 1 } as const;

const STEPS = [
  {
    label: 'Install',
    command: 'pip install "refgenie>=1.0.0a1"',
    note: 'Needs Python 3.11 or newer.',
  },
  {
    label: 'Pull an asset',
    command: 'refgenie pull hg38/fasta',
    note: 'Downloads it from the public server.',
  },
  {
    label: 'Use it',
    command: 'refgenie seek hg38/fasta',
    note: 'Prints the local path for your pipeline.',
  },
] as const;

function Command({ value, label }: { value: string; label?: string }) {
  return (
    <span className="rg-command flex items-center gap-2">
      <code className="rg-code rg-code--inline flex-1 overflow-x-auto">{value}</code>
      <CopyButton value={value} label={label ?? `Copy: ${value}`} />
    </span>
  );
}

/** One quiet line of counts, shown once they have all loaded. */
function Counts() {
  const genomes = useGenomes(COUNT_ONLY).data?.pagination.total;
  const groups = useAssetGroups(COUNT_ONLY).data?.pagination.total;
  const classes = useAssetClasses(COUNT_ONLY).data?.pagination.total;
  const recipes = useRecipes(COUNT_ONLY).data?.pagination.total;
  if ([genomes, groups, classes, recipes].some((n) => n === undefined)) return null;
  const n = (value: number | undefined) => (value ?? 0).toLocaleString();
  return (
    <p className="rg-muted text-sm">
      {n(genomes)} genomes · {n(groups)} asset groups · {n(classes)} asset classes ·{' '}
      {n(recipes)} recipes
    </p>
  );
}

export function LandingPage() {
  useDocumentTitle(undefined);
  const config = useUiConfig();
  const client = useApiClient();
  const navigate = useNavigate();
  const [query, setQuery] = useState('');

  const submit = (event: FormEvent) => {
    event.preventDefault();
    const q = query.trim();
    if (q !== '') navigate(`/genomes?q=${encodeURIComponent(q)}`);
  };

  // Root-relative for sub-path deployments; the absolute URL is for pasting into an assistant.
  const skillPath = `${config.root_path}/SKILL.md`;
  const skillUrl = new URL(skillPath, window.location.origin).href;
  const skillPrompt = `Read ${skillUrl}, then tell me which genomes this refgenie server has.`;

  // A public server offers MCP over HTTP; a local refgenie runs it over stdio.
  const isLocal = config.mode === 'local';
  const mcpCommand = isLocal
    ? 'claude mcp add refgenie refgenie-mcp'
    : `claude mcp add --transport http refgenie ${absoluteServerRoot(client.baseUrl)}/mcp`;

  return (
    <div className="flex flex-col gap-12">
      <section className="rg-hero text-center flex flex-col items-center gap-4">
        <img className="rg-hero__logo" src="refgenie_logo.svg" alt="" />
        <h1 className="text-4xl font-bold">Reference genome data management</h1>
        <p className="rg-hero__lede rg-muted text-lg">
          Refgenie downloads, builds and organizes the files your pipelines need for a
          reference genome: FASTA files, aligner indexes, annotations. Every machine
          gets the same files under the same names.
        </p>

        <form className="rg-hero__search flex gap-2 w-full" onSubmit={submit}>
          <label className="sr-only" htmlFor="landing-search">
            Search genomes
          </label>
          <input
            id="landing-search"
            className="rg-field__input"
            type="search"
            placeholder="Search genomes, e.g. hg38"
            value={query}
            onChange={(event) => setQuery(event.target.value)}
          />
          <button className="rg-btn rg-btn--primary" type="submit">
            Search
          </button>
        </form>

        <Counts />
      </section>

      <section className="flex flex-col gap-4">
        <h2 className="text-xl font-semibold">Get started</h2>
        <ol className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {STEPS.map((step, index) => (
            <li className="rg-step flex flex-col gap-2" key={step.command}>
              <span className="flex items-center gap-2">
                <span className="rg-step__index">{index + 1}</span>
                <span className="font-semibold">{step.label}</span>
              </span>
              <Command value={step.command} />
              <span className="rg-muted text-sm">{step.note}</span>
            </li>
          ))}
        </ol>
        <p className="rg-muted text-sm">
          The <ExternalLink href={PROJECT_LINKS.docs}>documentation</ExternalLink>{' '}
          covers building your own assets and running your own server.
        </p>
      </section>

      <section className="flex flex-col gap-4">
        <h2 className="text-xl font-semibold">Use it with an AI assistant</h2>

        <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
          <div className="rg-card rg-card__body flex flex-col gap-2">
            <span className="font-semibold">Connect the MCP server</span>
            <span className="rg-muted text-sm">
              {isLocal
                ? 'Lets an assistant look up your local genomes and assets. It cannot change anything.'
                : 'Lets an assistant look up the genomes and assets on this server. It cannot change anything.'}
            </span>
            <Command value={mcpCommand} />
          </div>

          <div className="rg-card rg-card__body flex flex-col gap-2">
            <span className="font-semibold">Or give it the API guide</span>
            <span className="rg-muted text-sm">
              {/* Plain <a>: /SKILL.md is a static file, not a client route. */}
              <a className="rg-link" href={skillPath}>
                SKILL.md
                <Icon name="code" className="rg-icon--trailing" />
              </a>{' '}
              explains this server&rsquo;s API in one page. Paste this into a chat:
            </span>
            <Command value={skillPrompt} label="Copy prompt" />
          </div>
        </div>

        <p className="rg-muted text-sm">
          To have an assistant run the command-line tool, give it the{' '}
          <ExternalLink href={PROJECT_LINKS.toolSkill}>CLI guide</ExternalLink>. More in{' '}
          <ExternalLink href={PROJECT_LINKS.useWithAi}>
            Use refgenie with AI
          </ExternalLink>
          .
        </p>
      </section>
    </div>
  );
}
