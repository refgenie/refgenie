import { ExternalLink } from '../common/ExternalLink';
import { useDocumentTitle } from '../../hooks/useDocumentTitle';
import { CollectionLogo } from './CollectionLogo';
import type { CollectionConfig } from '../../types/ui';

/** Head of a collection page: logo slot, name, blurb, external links and the note. */
export function CollectionHero({ collection }: { collection: CollectionConfig }) {
  useDocumentTitle(collection.name);
  return (
    <header className="rg-minihero rg-collection-hero">
      <div className="flex items-center gap-4">
        <CollectionLogo collection={collection} size="lg" />
        <h1 className="rg-minihero__title">{collection.name}</h1>
      </div>
      <p className="rg-minihero__lede">{collection.blurb}</p>
      {collection.links.length > 0 && (
        <ul className="rg-collection-hero__links" aria-label="External links">
          {collection.links.map((link) => (
            <li key={link.url}>
              <ExternalLink href={link.url}>{link.label}</ExternalLink>
            </li>
          ))}
        </ul>
      )}
      {collection.note && <p className="rg-muted text-xs">{collection.note}</p>}
    </header>
  );
}
