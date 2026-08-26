import { useUiConfig } from '../hooks/useUiConfig';
import { useApiClient } from '../hooks/useApiClient';
import { useConfigurations } from '../hooks/queries/useConfigurations';
import { DescriptionList } from '../components/common/DescriptionList';
import { ExternalLink } from '../components/common/ExternalLink';
import { CAPABILITY_KEYS } from '../types/ui';

export function AboutPage() {
  const config = useUiConfig();
  const client = useApiClient();
  // The configuration panel is a local-mode surface: the public server pages
  // never had it, and it names local filesystem paths.
  const configurations = useConfigurations({}, { enabled: config.mode === 'local' });
  const current = configurations.data?.items?.[0];

  return (
    <div className="flex flex-col gap-8">
      <h1 className="text-3xl font-bold">About</h1>

      <section>
        <h2 className="text-xl font-semibold mb-4">This instance</h2>
        <DescriptionList
          items={[
            { term: 'Service', value: config.service_name },
            { term: 'Mode', value: config.mode },
            { term: 'refgenie version', value: config.refgenie_version },
            { term: 'API base', value: <code className="rg-code rg-code--inline">{client.baseUrl}</code> },
            {
              term: 'Web UI build',
              value: config.web_ui?.commit
                ? `${config.web_ui.commit}${config.web_ui.dirty ? ' (dirty)' : ''}`
                : '',
            },
          ]}
        />
      </section>

      <section>
        <h2 className="text-xl font-semibold mb-4">Capabilities</h2>
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
      </section>

      {current && (
        <section>
          <h2 className="text-xl font-semibold mb-4">Configuration</h2>
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
        </section>
      )}

      <section>
        <h2 className="text-xl font-semibold mb-4">API</h2>
        <ul className="flex flex-col gap-1 text-sm">
          <li>
            <ExternalLink href={config.links.openapi ?? '/openapi.json'}>
              OpenAPI schema
            </ExternalLink>
          </li>
          <li>
            <ExternalLink href={`${config.root_path}/docs`}>Interactive API docs</ExternalLink>
          </li>
          {config.links.docs && (
            <li>
              <ExternalLink href={config.links.docs}>Documentation</ExternalLink>
            </li>
          )}
          {config.links.github && (
            <li>
              <ExternalLink href={config.links.github}>Source code</ExternalLink>
            </li>
          )}
        </ul>
      </section>
    </div>
  );
}
