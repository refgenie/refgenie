import { fastaDownloadEnabled } from '../services/config';
import type { FhrMetadata, GenomeDetailResponse } from '../types/api';
import { defaultRegion } from './parseRegion';

/** True if there is FHR metadata to show. */
export function hasFhr(fhr: FhrMetadata | null | undefined): fhr is FhrMetadata {
  return !!fhr && Object.keys(fhr).length > 0;
}

/** True if the genome page can offer a FASTA download for this genome. */
export function fastaDownloadAvailable(genome: GenomeDetailResponse): boolean {
  return fastaDownloadEnabled() && !!genome.refget_store_url;
}

export const GENOME_TABS = [
  'sequences',
  'extract',
  'download',
  'assets',
  'compare',
  'details',
] as const;
export type GenomeTab = (typeof GENOME_TABS)[number];

export function isGenomeTab(value: string | null): value is GenomeTab {
  return (GENOME_TABS as readonly string[]).includes(value ?? '');
}

/**
 * The tab to open when the URL names none: Sequences if the genome has
 * sequence data, else Assets if it has any, else Details. `null` while the
 * answer still depends on a request in flight.
 */
export function defaultGenomeTab({
  assetCount,
  hasSequences,
}: {
  assetCount: number | null;
  hasSequences: boolean | null;
}): GenomeTab | null {
  if (hasSequences === null) return null;
  if (hasSequences) return 'sequences';
  if (assetCount === null) return null;
  return assetCount > 0 ? 'assets' : 'details';
}

/** The Extract tab's link for a sequence: its first 1,000 bases. */
export function extractHref(row: { name: string; length: number }): string {
  const params = new URLSearchParams({ tab: 'extract', region: defaultRegion(row) });
  return `?${params.toString()}`;
}
