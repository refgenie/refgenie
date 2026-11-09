import { DescriptionList } from '../common/DescriptionList';
import { ExternalLink } from '../common/ExternalLink';
import { JsonBlock } from '../common/JsonBlock';
import type { DescriptionItem } from '../common/DescriptionList';
import type { FhrMetadata } from '../../types/api';
import { hasFhr } from '../../utils/genomePage';

export interface FhrPanelProps {
  fhr: FhrMetadata | null | undefined;
  /** Plain heading and content instead of a collapsed card. */
  flat?: boolean;
}

function names(people: FhrMetadata['metadataAuthor']): string | undefined {
  if (!people || people.length === 0) return undefined;
  const joined = people
    .map((person) => person?.name)
    .filter((name): name is string => !!name)
    .join(', ');
  return joined || undefined;
}

/**
 * The FHR sidecar with a raw-JSON toggle: a collapsed card, or with `flat`,
 * a plain section.
 */
export function FhrPanel({ fhr, flat = false }: FhrPanelProps) {
  if (!hasFhr(fhr)) return null;

  const stats = fhr.vitalStats ?? {};
  const statItems: DescriptionItem[] = (
    [
      ['N50', stats.N50],
      ['L50', stats.L50],
      ['L90', stats.L90],
      ['Total base pairs', stats.totalBasePairs],
      ['Contigs', stats.numberContigs],
      ['Scaffolds', stats.numberScaffolds],
      ['Read technology', stats.readTechnology],
    ] as const
  )
    .filter(([, value]) => value !== null && value !== undefined)
    .map(([term, value]) => ({ term, value: String(value) }));

  const provenance: DescriptionItem[] = [
    { term: 'Date created', value: fhr.dateCreated ?? '' },
    { term: 'Metadata authors', value: names(fhr.metadataAuthor) ?? '' },
    { term: 'Assembly authors', value: names(fhr.assemblyAuthor) ?? '' },
    { term: 'License', value: fhr.license ?? '' },
    {
      term: 'Publication',
      value: fhr.scholarlyArticle ? (
        <ExternalLink href={fhr.scholarlyArticle}>{fhr.scholarlyArticle}</ExternalLink>
      ) : (
        ''
      ),
    },
    { term: 'Funding', value: fhr.funding ?? '' },
    {
      term: 'Related links',
      value:
        fhr.relatedLink && fhr.relatedLink.length > 0 ? (
          <span className="flex flex-col gap-1">
            {fhr.relatedLink.map((link) => (
              <ExternalLink key={link} href={link}>
                {link}
              </ExternalLink>
            ))}
          </span>
        ) : (
          ''
        ),
    },
  ];

  const body = (
    <>
      {statItems.length > 0 && (
        <section>
          <h3 className="text-sm font-semibold mb-2">Assembly stats</h3>
          <DescriptionList items={statItems} hideEmpty />
        </section>
      )}
      <section>
        <h3 className="text-sm font-semibold mb-2">Provenance and publication</h3>
        <DescriptionList items={provenance} hideEmpty />
      </section>
      <JsonBlock value={fhr} label="Raw metadata JSON" />
    </>
  );

  if (flat) {
    return (
      <section className="rg-flat-section flex flex-col gap-6">
        <h2>Metadata</h2>
        {body}
      </section>
    );
  }

  return (
    <details className="rg-flat-details">
      <summary className="rg-flat-details__summary">Metadata</summary>
      <div className="flex flex-col gap-6">{body}</div>
    </details>
  );
}
