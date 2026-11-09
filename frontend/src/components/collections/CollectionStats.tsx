import type { CollectionStats as Stats } from '../../utils/collections';
import { HAPLOTYPES } from '../../utils/haplotype';
import type { CollectionConfig } from '../../types/ui';

export interface CollectionStatsProps {
  collection: CollectionConfig;
  /** `undefined` while the genomes load. */
  stats: Stats | undefined;
}

/** Headline counts, like the GenomesPage summary cards. */
export function CollectionStats({ collection, stats }: CollectionStatsProps) {
  const cards: Array<[string, number | undefined]> = [
    ['Genomes', stats?.genomes],
    ['Species', stats?.species],
    ['Genomes with built assets', stats?.withAssets],
    ['Asset groups', stats?.assetGroups],
  ];
  const showHaplotypes = collection.columns.includes('haplotype') && stats;
  return (
    <div className="flex flex-col gap-2">
      <dl className="grid grid-cols-2 sm:grid-cols-4 gap-4">
        {cards.map(([term, value]) => (
          <div className="rg-card p-4" key={term}>
            <dt className="rg-kv__term">{term}</dt>
            <dd className="text-2xl font-semibold">{value ?? '—'}</dd>
          </div>
        ))}
      </dl>
      {showHaplotypes && (
        <p className="rg-muted text-sm" aria-label="Haplotype counts">
          {[...HAPLOTYPES, 'unspecified' as const]
            .filter((h) => stats.haplotypes[h] > 0)
            .map((h) => `${stats.haplotypes[h]} ${h}`)
            .join(' · ')}
        </p>
      )}
    </div>
  );
}
