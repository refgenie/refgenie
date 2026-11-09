import { Copyable } from '../common/Copyable';
import { DescriptionList } from '../common/DescriptionList';

export interface CollectionDigestsProps {
  /** Level-1 attribute digests: names, lengths, sequences, ... */
  level1: Record<string, string>;
}

/** The Details tab's per-attribute digests of the genome's sequence collection. */
export function CollectionDigests({ level1 }: CollectionDigestsProps) {
  const items = Object.entries(level1).map(([attribute, digest]) => ({
    term: attribute,
    value: (
      <Copyable value={digest} label={`Copy ${attribute} digest`} mono>
        {digest}
      </Copyable>
    ),
  }));

  return (
    <section className="rg-flat-section">
      <h2>Collection digests</h2>
      <DescriptionList items={items} />
    </section>
  );
}
