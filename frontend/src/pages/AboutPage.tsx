/**
 * `/about` — two concerns, one page.
 *
 * 1. Refgenie the project: what it is, where the docs are, who to cite, how to
 *    report a bug. Identical in both modes on purpose — the project is the same
 *    thing whether you are reading this on a public server or on your laptop.
 * 2. This instance: version, capabilities, configuration, endpoints. The only
 *    `mode` branch on the page is the configuration panel, which names local
 *    filesystem paths; everything else is capability-gated.
 *
 * Project URLs are constants from services/projectLinks.ts, NOT the `links`
 * block of /service-info, which the production server fills with the retired
 * readthedocs URL. Instance URLs go through instanceUrl() so they reach the API
 * origin on a cross-origin static deployment rather than the SPA's own.
 *
 * The landing page owns the "what does refgenie do for me" pitch and the three
 * install commands. Do not repeat either here.
 */

import { Link } from 'react-router-dom';
import { useUiConfig } from '../hooks/useUiConfig';
import { useApiClient } from '../hooks/useApiClient';
import { useConfigurations } from '../hooks/queries/useConfigurations';
import { MiniHero } from '../components/layout/MiniHero';
import { DescriptionList } from '../components/common/DescriptionList';
import { ApiLink } from '../components/common/ApiLink';
import { ExternalLink } from '../components/common/ExternalLink';
import { CopyButton } from '../components/common/CopyButton';
import { Icon } from '../components/common/Icon';
import { instanceUrl } from '../services/instanceLinks';
import { CITATIONS, DOC_GROUPS, ECOSYSTEM, PROJECT_LINKS } from '../services/projectLinks';
import type { LinkEntry } from '../services/projectLinks';
import { formatTimestamp } from '../utils/time';
import { CAPABILITY_KEYS } from '../types/ui';

/**
 * An endpoint on this instance.
 *
 * `kind` picks the link component, and with it the glyph. `json` is a raw
 * document this service serves — `ApiLink`'s `</>`. `page` is an HTML API
 * explorer, which is a site you go to rather than a payload you read, so it
 * keeps the outward arrow. That is the same split the OpenAPI schema and the
 * interactive docs already had on this page before it was rewritten.
 */
interface Endpoint {
  label: string;
  path: string;
  kind: 'json' | 'page';
}

function LinkCard({ entry }: { entry: LinkEntry }) {
  return (
    <ExternalLink
      className="rg-card rg-card__body rg-link--plain no-underline flex flex-col gap-2"
      href={entry.href}
    >
      <span className="font-semibold">{entry.label}</span>
      <span className="rg-muted text-sm">{entry.blurb}</span>
    </ExternalLink>
  );
}

export function AboutPage() {
  const config = useUiConfig();
  const client = useApiClient();
  // The configuration panel is a local-mode surface: the public server pages
  // never had it, and it names local filesystem paths.
  const configurations = useConfigurations({}, { enabled: config.mode === 'local' });
  const current = configurations.data?.items?.[0];

  const endpoints: Endpoint[] = [
    { label: 'OpenAPI schema', path: '/openapi.json', kind: 'json' },
    { label: 'Interactive API docs', path: '/docs', kind: 'page' },
    { label: 'API reference (ReDoc)', path: '/redoc', kind: 'page' },
    { label: 'Service info document', path: '/service-info', kind: 'json' },
  ];
  // Capability-gated, never mode-gated: a local dash mounts neither router, and
  // reports both flags false, so these rows vanish on their own.
  if (config.capabilities.seqcol) {
    endpoints.push({
      label: 'Sequence collections service info',
      path: '/seqcol/service-info',
      kind: 'json',
    });
  }
  if (config.capabilities.drs) {
    // Mounted at /ga4gh/drs (and /v4/ga4gh/drs) in refgenie/server/main.py.
    // There is no /v1 segment: /ga4gh/drs/v1/service-info is a 404.
    endpoints.push({
      label: 'GA4GH DRS service info',
      path: '/ga4gh/drs/service-info',
      kind: 'json',
    });
  }

  return (
    <div className="flex flex-col gap-12">
      <MiniHero
        title="About refgenie"
        documentTitle="About"
        lede={
          <>
            Refgenie the open-source project — its documentation, its papers, its source — and
            the particular instance serving you this page, with the version it runs and the
            operations it will let you perform.
          </>
        }
        // The jump nav sits on the title's line rather than under the lede: it
        // is navigation into this page, which is what the actions slot is for.
        actions={
          <nav className="flex flex-wrap gap-4 text-sm" aria-label="On this page">
            <a className="rg-link" href="#project">
              The project
            </a>
            <a className="rg-link" href="#instance">
              This instance
            </a>
          </nav>
        }
      />

      {/* ============================ THE PROJECT ============================ */}
      <section id="project" aria-labelledby="project-heading" className="flex flex-col gap-8">
        <div className="flex flex-col gap-3">
          <h2 id="project-heading" className="text-2xl font-semibold">
            Refgenie, the project
          </h2>
          <p className="rg-muted">
            Refgenie is open-source software for managing reference genome assets. It is a
            Python package with a command-line interface, a REST API, and this web interface,
            all from one install. Genomes are identified by GA4GH sequence-collection digests
            derived from the sequences themselves, so two people can tell whether they are
            working from the same assembly without comparing files. It is developed in the open
            by the Sheffield Lab at the University of Virginia and released under the BSD
            2-Clause license.
          </p>
          <p className="rg-muted text-sm">
            New here? The{' '}
            <Link className="rg-link" to="/">
              home page
            </Link>{' '}
            has the three commands that get you a genome asset on disk.
          </p>
        </div>

        <div className="flex flex-col gap-6">
          <h3 className="text-xl font-semibold">Documentation</h3>
          {DOC_GROUPS.map((group) => (
            <div className="flex flex-col gap-3" key={group.heading}>
              <h4 className="font-semibold text-sm rg-muted">{group.heading}</h4>
              <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
                {group.items.map((entry) => (
                  <LinkCard entry={entry} key={entry.href} />
                ))}
              </div>
            </div>
          ))}
        </div>

        <div className="flex flex-col gap-4">
          <h3 className="text-xl font-semibold">Source, issues, and license</h3>
          <DescriptionList
            items={[
              {
                term: 'Source code',
                value: (
                  <ExternalLink href={PROJECT_LINKS.github}>
                    github.com/refgenie/refgenie
                  </ExternalLink>
                ),
              },
              {
                term: 'Report a problem',
                value: (
                  <ExternalLink href={PROJECT_LINKS.issues}>
                    Open an issue on GitHub
                  </ExternalLink>
                ),
              },
              {
                term: 'Package',
                value: (
                  <ExternalLink href={PROJECT_LINKS.pypi}>refgenie on PyPI</ExternalLink>
                ),
              },
              {
                term: 'License',
                value: <ExternalLink href={PROJECT_LINKS.license}>BSD 2-Clause</ExternalLink>,
              },
              {
                term: 'Data channel registry',
                value: (
                  <ExternalLink href={PROJECT_LINKS.registry}>refgenie-registry</ExternalLink>
                ),
              },
              {
                term: 'Documentation source',
                value: (
                  <ExternalLink href={PROJECT_LINKS.docsSource}>refgenie-docs</ExternalLink>
                ),
              },
            ]}
          />
        </div>

        <div className="flex flex-col gap-4">
          <h3 className="text-xl font-semibold">How to cite refgenie</h3>
          <p className="rg-muted text-sm">
            If refgenie helped with work you are publishing, please cite it.
          </p>
          <ol className="flex flex-col gap-4">
            {CITATIONS.map((citation) => (
              <li className="rg-card rg-card__body flex flex-col gap-2" key={citation.doi}>
                <span>
                  {citation.authors} ({citation.year}).{' '}
                  <ExternalLink href={citation.doi}>{citation.title}</ExternalLink>.{' '}
                  <span className="rg-muted">{citation.venue}</span>
                </span>
                <span className="rg-muted text-sm">{citation.note}</span>
                <span className="flex items-center gap-2">
                  <code className="rg-code rg-code--inline flex-1 overflow-x-auto">
                    {citation.doi}
                  </code>
                  <CopyButton
                    value={citation.plain}
                    label={`Copy citation: ${citation.title}`}
                  />
                </span>
              </li>
            ))}
          </ol>
        </div>

        <div className="flex flex-col gap-4">
          <h3 className="text-xl font-semibold">The wider ecosystem</h3>
          <p className="rg-muted text-sm">
            Refgenie sits on top of the GA4GH standards for identifying sequences and
            assemblies, and shares its identity layer with the refget project.
          </p>
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
            {ECOSYSTEM.map((entry) => (
              <LinkCard entry={entry} key={entry.href} />
            ))}
          </div>
        </div>
      </section>

      {/* =========================== THIS INSTANCE =========================== */}
      <section id="instance" aria-labelledby="instance-heading" className="flex flex-col gap-8">
        <div className="flex flex-col gap-3">
          <h2 id="instance-heading" className="text-2xl font-semibold">
            This instance
          </h2>
          <p className="rg-muted">
            Everything below describes the {config.service_name} you are connected to right
            now, not refgenie in general. A different instance holds different genomes and
            allows different operations.
          </p>
          {config.degraded && (
            <p className="rg-muted text-sm">
              This page could not read{' '}
              <code className="rg-code rg-code--inline">/service-info</code>, so the facts
              below are defaults rather than what the server actually reports.
            </p>
          )}
        </div>

        <div className="flex flex-col gap-4">
          <h3 className="text-xl font-semibold">Identity</h3>
          <DescriptionList
            items={[
              { term: 'Service', value: config.service_name },
              { term: 'Mode', value: config.mode },
              { term: 'refgenie version', value: config.refgenie_version },
              {
                term: 'API base',
                value: <code className="rg-code rg-code--inline">{client.baseUrl}</code>,
              },
              {
                term: 'Web UI build',
                value: config.web_ui?.commit
                  ? `${config.web_ui.commit}${config.web_ui.dirty ? ' (dirty)' : ''}`
                  : '',
              },
              {
                term: 'Web UI built',
                value: config.web_ui?.built_at ? formatTimestamp(config.web_ui.built_at) : '',
              },
            ]}
          />
        </div>

        <div className="flex flex-col gap-4">
          <h3 className="text-xl font-semibold">Capabilities</h3>
          <p className="rg-muted text-sm">
            Every action in this interface is gated on one of these flags. A public server
            turns the command flags off; a local dashboard turns them on.
          </p>
          <div className="rg-table__wrap">
            <table className="rg-table">
              <caption className="sr-only">Capabilities of this instance</caption>
              <thead className="rg-table__head">
                <tr>
                  <th scope="col" className="rg-table__cell">
                    Capability
                  </th>
                  <th scope="col" className="rg-table__cell">
                    Enabled
                  </th>
                </tr>
              </thead>
              <tbody>
                {CAPABILITY_KEYS.map((key) => (
                  <tr className="rg-table__row" key={key}>
                    <td className="rg-table__cell">
                      <code className="rg-code rg-code--inline">{key}</code>
                    </td>
                    <td className="rg-table__cell">{config.capabilities[key] ? 'yes' : 'no'}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>

        {current && (
          <div className="flex flex-col gap-4">
            <h3 className="text-xl font-semibold">Configuration</h3>
            <DescriptionList
              items={[
                { term: 'Config version', value: String(current.version) },
                { term: 'Genome folder', value: current.genome_folder },
                { term: 'Genome stage folder', value: current.genome_stage_folder ?? '' },
                {
                  term: 'Servers',
                  value: current.servers.length ? (
                    <ul className="flex flex-col gap-1">
                      {current.servers.map((server) => (
                        <li key={server}>
                          <code className="rg-code rg-code--inline">{server}</code>
                        </li>
                      ))}
                    </ul>
                  ) : (
                    ''
                  ),
                },
              ]}
            />
          </div>
        )}

        <div className="flex flex-col gap-4">
          <h3 className="text-xl font-semibold">Endpoints</h3>
          <p className="rg-muted text-sm">
            Everything this interface shows comes from these; you can call them directly.
          </p>
          <ul className="rg-list flex flex-col gap-1">
            <li>
              {/* A plain <a>, and NOT an ApiLink: /SKILL.md is a static file on
                  this origin, not a client route, so the SPA router must never
                  see it, and it deliberately stays in this tab. It still earns
                  the `code` glyph — it is a raw document this instance serves,
                  the same class as the OpenAPI schema below. The visible label
                  already names the file, so no sr-only text is needed. */}
              <a className="rg-link" href={`${config.root_path}/SKILL.md`}>
                Agent capability doc (SKILL.md)
                <Icon name="code" className="rg-icon--trailing" />
              </a>
            </li>
            {endpoints.map((endpoint) => {
              const Anchor = endpoint.kind === 'json' ? ApiLink : ExternalLink;
              return (
                <li key={endpoint.path}>
                  <Anchor href={instanceUrl(client.baseUrl, endpoint.path)}>
                    {endpoint.label}
                  </Anchor>{' '}
                  <code className="rg-code rg-code--inline">{endpoint.path}</code>
                </li>
              );
            })}
          </ul>
        </div>
      </section>
    </div>
  );
}
