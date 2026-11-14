import type { CollectionConfig } from '../../types/ui';
import { cn } from '../../utils/cn';
import { safeHttpUrl } from '../../utils/links';

/** Up to three initials from the name's capitalized words, e.g. "VGPA". */
function initials(name: string): string {
  const letters = name
    .split(/\s+/)
    .filter((word) => /^[A-Z]/.test(word))
    .map((word) => word[0]);
  return (letters.length ? letters : [name[0] ?? '?'])
    .slice(0, 3)
    .join('')
    .toUpperCase();
}

export interface CollectionLogoProps {
  collection: Pick<CollectionConfig, 'name' | 'short_name' | 'logo_url'>;
  size?: 'md' | 'lg';
}

/**
 * The collection's logo, or a neutral initials badge on the accent color.
 * No third-party logos are bundled; a deployment supplies `logo_url`.
 */
export function CollectionLogo({ collection, size = 'md' }: CollectionLogoProps) {
  const className = cn('rg-collection-logo', size === 'lg' && 'rg-collection-logo--lg');
  const logo = safeHttpUrl(collection.logo_url);
  if (logo) {
    return (
      <img
        className={className}
        src={logo}
        alt={`${collection.name} logo`}
      />
    );
  }
  return (
    <span className={cn(className, 'rg-collection-logo--badge')} aria-hidden="true">
      {collection.short_name ?? initials(collection.name)}
    </span>
  );
}
