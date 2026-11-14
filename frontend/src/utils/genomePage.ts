import { fastaDownloadEnabled } from '../services/config';
import type { AssetResponse, FhrMetadata, GenomeDetailResponse } from '../types/api';
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

/** The Extract tab's link for a sequence: its first 1,000 bases. */
export function extractHref(row: { name: string; length: number }): string {
  const params = new URLSearchParams({ tab: 'extract', region: defaultRegion(row) });
  return `?${params.toString()}`;
}

/**
 * The size of what `refgenie pull <alias>/fasta` fetches: the `fasta` group's
 * default asset, else its first. Null when there is none or its size is unknown.
 */
export function fastaPullSize(assets: readonly AssetResponse[] | undefined): number | null {
  const group = assets?.filter((asset) => asset.asset_group_name === 'fasta') ?? [];
  const asset = group.find((a) => a.is_default) ?? group[0];
  return asset?.size ?? null;
}
