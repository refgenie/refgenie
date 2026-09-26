import { Icon } from './Icon';
import type { IconName } from './Icon';

export interface ExternalLinkProps {
  href: string;
  children: React.ReactNode;
  className?: string;
  /**
   * The affordance glyph. `external` (the default) means the link leaves this
   * site. `code` means it opens this instance's own raw JSON — same tab
   * behaviour, different promise, so a different glyph. `none` is for the rare
   * place the indicator would be pure decoration.
   */
  icon?: IconName | 'none';
  /**
   * Screen-reader suffix. The icon is `aria-hidden`, so the fact it carries has
   * to reach the accessibility tree as text. Pass `''` to suppress.
   */
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
