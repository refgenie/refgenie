import { useUiConfig } from '../../hooks/useUiConfig';
import { ExternalLink } from '../common/ExternalLink';

export function Footer() {
  const config = useUiConfig();
  const links = Object.entries(config.links).filter(([, href]) => !!href) as Array<
    [string, string]
  >;

  return (
    <footer className="rg-layout__footer p-4 flex flex-wrap items-center justify-between gap-4">
      <span>
        {config.service_name} · refgenie {config.refgenie_version}
      </span>
      <span className="flex flex-wrap gap-4">
        {links.map(([name, href]) => (
          <ExternalLink key={name} href={href}>
            {name}
          </ExternalLink>
        ))}
      </span>
    </footer>
  );
}
