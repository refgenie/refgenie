/**
 * Site footer. Project links are frontend constants rather than the server's
 * `config.links`; the instance link is built from the API base.
 */

import { Link } from 'react-router-dom';
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
        <Link className="rg-link" to="/about">
          About
        </Link>
        <ExternalLink href={PROJECT_LINKS.github}>Source</ExternalLink>
        <ExternalLink href={PROJECT_LINKS.issues}>Report a problem</ExternalLink>
        <ApiLink href={instanceUrl(client.baseUrl, '/openapi.json')}>
          API schema
        </ApiLink>
      </span>
    </footer>
  );
}
