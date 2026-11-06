import { Icon } from './Icon';
import type { IconName } from './Icon';

export interface ExternalLinkProps {
  href: string;
  children: React.ReactNode;
  className?: string;
  /** Glyph: `external` (default) leaves the site, `code` opens this instance's raw JSON, `none` hides it. */
  icon?: IconName | 'none';
  /** Screen-reader suffix standing in for the hidden icon. Pass `''` to suppress. */
  hint?: string;
  title?: string;
}

export function ExternalLink({
  href,
  children,
  className = 'rg-link',
  icon = 'external',
  hint = 'opens in a new tab',
  title,
}: ExternalLinkProps) {
  return (
    <a
      href={href}
      className={className}
      target="_blank"
      rel="noopener noreferrer"
      title={title}
    >
      {children}
      {icon !== 'none' && <Icon name={icon} className="rg-icon--trailing" />}
      {hint !== '' && <span className="sr-only"> ({hint})</span>}
    </a>
  );
}
