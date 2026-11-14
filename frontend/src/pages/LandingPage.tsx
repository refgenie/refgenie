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
import { Copyable } from '../components/common/Copyable';
import { ExternalLink } from '../components/common/ExternalLink';
import { Icon } from '../components/common/Icon';
import { Badge } from '../components/common/Badge';
import { CopyButton } from '../components/common/CopyButton';
import { installCommand } from '../utils/install';
import { PROJECT_LINKS } from '../services/projectLinks';
import { absoluteServerRoot } from '../services/instanceLinks';

/** Count-only request: one row per page. */
const COUNT_ONLY = { pageSize: 1 } as const;

/**
 * The quick start after install, the same everywhere: a first pull that is
 * small (hg38's FASTA index and chromosome sizes), not the 3 GB FASTA.
 */
const LATER_STEPS = [
  {
    label: 'Set up',
    command: 'refgenie init',
    note: 'Creates your refgenie folder and connects it to the public server.',
  },
  {
    label: 'Pull an asset',
    command: 'refgenie pull hg38/fasta_index',
    note: 'About 12 KB: the FASTA index (.fai) and chromosome sizes for hg38.',
  },
  {
    label: 'Use it',
    command: 'refgenie seek hg38/fasta_index.chrom_sizes',
    note: 'Prints the local path of the file, for your pipeline.',
  },
] as const;

function Command({ value, label }: { value: string; label?: string }) {
  return (
    <Copyable value={value} label={label ?? `Copy: ${value}`} className="rg-command">
      <code className="rg-code rg-code--inline">{value}</code>
    </Copyable>
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
        <p className="rg-hero__note text-sm">
          <Badge>Alpha</Badge> Refgenie 1.0 is an alpha release: it works, but commands
          and data can still change, so pin the version in production pipelines.{' '}
          <ExternalLink href={PROJECT_LINKS.refgenieDocs}>Details</ExternalLink>
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
        <ol className="grid grid-cols-1 md:grid-cols-2 gap-6">
          {[
            {
              label: 'Install',
              command: installCommand(config.refgenie_version),
              note: 'Needs Python 3.11 or newer.',
            },
            ...LATER_STEPS,
          ].map((step, index) => (
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
          Then pull what your pipeline needs, such as{' '}
          <code className="rg-code rg-code--inline">hg38/bowtie2_index</code>. The{' '}
          <ExternalLink href={PROJECT_LINKS.docs}>documentation</ExternalLink> covers
          building your own assets and running your own server.
        </p>
      </section>

      {/* Kept short: most visitors come for the genomes and the CLI. */}
      <section className="flex flex-col gap-2">
        <h2 className="text-lg font-semibold">Use it with an AI assistant</h2>
        <p className="rg-muted text-sm">
          {isLocal
            ? 'Connect the MCP server to let an assistant look up your local genomes and assets. It cannot change anything.'
            : 'Connect the MCP server to let an assistant look up the genomes and assets on this server. It cannot change anything.'}
        </p>
        <div className="rg-landing__ai-command">
          <Command value={mcpCommand} />
        </div>
        <p className="rg-muted text-sm">
          Or give it{' '}
          {/* Plain <a>: /SKILL.md is a static file, not a client route. */}
          <a className="rg-link" href={skillPath}>
            SKILL.md
            <Icon name="code" className="rg-icon--trailing" />
          </a>
          , this server&rsquo;s API in one page{' '}
          <span className="rg-landing__prompt">
            (copy a starter prompt <CopyButton value={skillPrompt} label="Copy prompt" />)
          </span>
          , or the <ExternalLink href={PROJECT_LINKS.toolSkill}>CLI guide</ExternalLink> to
          run the command-line tool. More in{' '}
          <ExternalLink href={PROJECT_LINKS.useWithAi}>Use refgenie with AI</ExternalLink>.
        </p>
      </section>
    </div>
  );
}
