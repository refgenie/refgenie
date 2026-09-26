import { ExternalLink } from './ExternalLink';

export interface ApiLinkProps {
  href: string;
  children: React.ReactNode;
  className?: string;
}

/**
 * A link to a raw API response served by this instance.
 *
 * The transport is `ExternalLink`'s — it opens a new tab and leaves the SPA —
 * but the destination is this service's own JSON, not somebody else's site.
 * Keeping the two apart is what keeps the outward arrow honest: an arrow should
 * mean "another site", and if it also means "our own /openapi.json" it means
 * nothing.
 */
export function ApiLink({ href, children, className }: ApiLinkProps) {
  return (
    <ExternalLink href={href} className={className} icon="code">
      {children}
    </ExternalLink>
  );
}
