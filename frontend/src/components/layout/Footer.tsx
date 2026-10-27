/**
 * Project URLs are frontend constants, not `config.links` from the server:
 * iterating that would put raw wire keys (`docs`, `github`, `openapi`) on
 * screen as link text and ship whatever the server says `docs` is, which in
 * production is the retired readthedocs site. The one instance link is derived
 * from the API base, so it survives a cross-origin deployment.
 */

import { useUiConfig } from '../../hooks/useUiConfig';
import { useApiClient } from '../../hooks/useApiClient';
import { ApiLink } from '../common/ApiLink';
import { ExternalLink } from '../common/ExternalLink';
import { instanceUrl } from '../../services/instanceLinks';
import { PROJECT_LINKS } from '../../services/projectLinks';

export function Footer() {
  const config = useUiConfig();
  const client = useApiClient();

  return (
    <footer className="rg-layout__footer p-4 flex flex-wrap items-center justify-between gap-4">
      <span>
        {config.service_name} · refgenie {config.refgenie_version}
      </span>
      <span className="flex flex-wrap gap-4">
        <ExternalLink href={PROJECT_LINKS.docs}>Documentation</ExternalLink>
        <ExternalLink href={PROJECT_LINKS.github}>Source</ExternalLink>
        <ExternalLink href={PROJECT_LINKS.issues}>Report a problem</ExternalLink>
        <ApiLink href={instanceUrl(client.baseUrl, '/openapi.json')}>API schema</ApiLink>
      </span>
    </footer>
  );
}
