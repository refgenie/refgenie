/** `/collections`: one card per collection the server declares. */

import { MiniHero } from '../components/layout/MiniHero';
import { ScopedLink } from '../components/common/ScopedLink';
import { EmptyState } from '../components/common/states';
import { CollectionLogo } from '../components/collections/CollectionLogo';
import { useAccent } from '../components/collections/useAccent';
import { useCollectionCount, useCollections } from '../hooks/queries/useCollections';
import type { CollectionConfig } from '../types/ui';

/** The blurb's first sentence, for a card. */
function firstSentence(text: string): string {
  const match = /^.*?[.!?](?=\s|$)/.exec(text.trim());
  return match ? match[0] : text;
}

function CollectionCard({ collection }: { collection: CollectionConfig }) {
  const count = useCollectionCount(collection);
  const accentRef = useAccent<HTMLLIElement>(collection.accent_color);
  return (
    <li className="rg-card rg-collection-card" ref={accentRef}>
      <div className="rg-card__body flex flex-col gap-3">
        <div className="flex items-center gap-3">
          <CollectionLogo collection={collection} />
          <h2 className="rg-card__title">
            <ScopedLink className="rg-link" to={`/collections/${collection.id}`}>
              {collection.name}
            </ScopedLink>
          </h2>
        </div>
        <p className="text-sm">{firstSentence(collection.blurb)}</p>
        <p className="rg-muted text-sm">
          {count.data === undefined ? '—' : count.data} genome
          {count.data === 1 ? '' : 's'}
        </p>
      </div>
    </li>
  );
}

export function CollectionsPage() {
  const collections = useCollections();
  return (
    <div className="flex flex-col gap-6">
      <MiniHero
        title="Collections"
        lede="Curated groups of genomes, each drawn from one or more refget stores on this server."
      />
      {collections.length === 0 ? (
        <EmptyState message="This server has no collections configured." />
      ) : (
        <ul className="rg-collection-grid grid grid-cols-1 md:grid-cols-2 gap-4">
          {collections.map((collection) => (
            <CollectionCard key={collection.id} collection={collection} />
          ))}
        </ul>
      )}
    </div>
  );
}
