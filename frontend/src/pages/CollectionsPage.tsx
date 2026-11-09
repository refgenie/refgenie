/** `/collections`: one row per collection the server declares. */

import { MiniHero } from '../components/layout/MiniHero';
import { ScopedLink } from '../components/common/ScopedLink';
import { EmptyState } from '../components/common/states';
import { CollectionLogo } from '../components/collections/CollectionLogo';
import { useAccent } from '../components/collections/useAccent';
import { useCollectionCount, useCollections } from '../hooks/queries/useCollections';
import type { CollectionConfig } from '../types/ui';

/** The blurb's first sentence, for a row. */
function firstSentence(text: string): string {
  const match = /^.*?[.!?](?=\s|$)/.exec(text.trim());
  return match ? match[0] : text;
}

function plural(n: number, word: string): string {
  return `${n} ${word}${n === 1 ? '' : 's'}`;
}

function CollectionRow({ collection }: { collection: CollectionConfig }) {
  const count = useCollectionCount(collection);
  const accentRef = useAccent<HTMLLIElement>(collection.accent_color);
  const genomes =
    count.data === undefined ? '… genomes' : plural(count.data, 'genome');
  return (
    <li className="rg-collection-row" ref={accentRef}>
      <CollectionLogo collection={collection} />
      <div className="flex flex-col gap-1 min-w-0">
        <h2 className="text-lg font-semibold">
          <ScopedLink className="rg-link" to={`/collections/${collection.id}`}>
            {collection.name}
          </ScopedLink>
        </h2>
        <p className="text-sm">{firstSentence(collection.blurb)}</p>
        <p className="rg-muted text-sm">
          {genomes} · {plural(collection.stores.length, 'store')}
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
        <ul className="rg-collection-list">
          {collections.map((collection) => (
            <CollectionRow key={collection.id} collection={collection} />
          ))}
        </ul>
      )}
    </div>
  );
}
