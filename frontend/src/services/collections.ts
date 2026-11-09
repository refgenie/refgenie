/** Coerce the wire `collections` value from `/service-info` into typed configs. */

import type {
  CollectionColumn,
  CollectionConfig,
  CollectionLink,
  GenomeLinkTemplate,
} from '../types/ui';

const COLUMNS: readonly CollectionColumn[] = [
  'assembly_accession',
  'haplotype',
  'assembly_level',
  'assembly_source',
];

const HEX_COLOR = /^#[0-9a-fA-F]{6}$/;

function isRecord(value: unknown): value is Record<string, unknown> {
  return typeof value === 'object' && value !== null && !Array.isArray(value);
}

function optionalString(value: unknown): string | null {
  return typeof value === 'string' && value !== '' ? value : null;
}

function links(value: unknown): CollectionLink[] {
  if (!Array.isArray(value)) return [];
  return value.filter(
    (l): l is CollectionLink =>
      isRecord(l) && typeof l.label === 'string' && typeof l.url === 'string',
  );
}

function genomeLinks(value: unknown): GenomeLinkTemplate[] {
  if (!Array.isArray(value)) return [];
  return value.filter(
    (l): l is GenomeLinkTemplate =>
      isRecord(l) && typeof l.label === 'string' && typeof l.template === 'string',
  );
}

function coerceOne(raw: unknown): CollectionConfig | null {
  if (!isRecord(raw)) return null;
  const { id, name, blurb, accent_color, stores } = raw;
  if (typeof id !== 'string' || typeof name !== 'string' || typeof blurb !== 'string') {
    return null;
  }
  if (typeof accent_color !== 'string' || !HEX_COLOR.test(accent_color)) return null;
  if (
    !Array.isArray(stores) ||
    stores.length === 0 ||
    !stores.every((s) => typeof s === 'string' && s !== '')
  ) {
    return null;
  }
  const columns = Array.isArray(raw.columns)
    ? raw.columns.filter((c): c is CollectionColumn =>
        COLUMNS.includes(c as CollectionColumn),
      )
    : [];
  return {
    id,
    name,
    short_name: optionalString(raw.short_name),
    blurb,
    note: optionalString(raw.note),
    accent_color,
    logo_url: optionalString(raw.logo_url),
    stores: stores as string[],
    links: links(raw.links),
    genome_links: genomeLinks(raw.genome_links),
    columns,
  };
}

/** Keep only well-formed collection entries; anything else becomes `[]`. */
export function coerceCollections(raw: unknown): CollectionConfig[] {
  if (!Array.isArray(raw)) return [];
  return raw.map(coerceOne).filter((c): c is CollectionConfig => c !== null);
}
