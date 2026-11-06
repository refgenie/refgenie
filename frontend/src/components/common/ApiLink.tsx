import { ExternalLink } from './ExternalLink';

export interface ApiLinkProps {
  href: string;
  children: React.ReactNode;
  className?: string;
}

/**
 * A link to a raw API response on this instance. Uses the `code` glyph, since
 * the outward arrow is reserved for other sites.
 */
export function ApiLink({ href, children, className }: ApiLinkProps) {
  return (
    <ExternalLink href={href} className={className} icon="code">
      {children}
    </ExternalLink>
  );
}
