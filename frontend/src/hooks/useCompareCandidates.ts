import { useMemo } from 'react';
import { useGenomes } from './queries/useGenomes';
import { rankCompareCandidates } from '../utils/compareCandidates';
import type { CompareCandidate } from '../utils/compareCandidates';
import type { GenomeDetailResponse } from '../types/api';

/**
 * The genomes `genome` can be compared with: the others of its species,
 * likeliest first. Asks the server for that species only (by taxon id, else
 * by species name). `null` while loading.
 */
export function useCompareCandidates(
  genome: GenomeDetailResponse | undefined,
): CompareCandidate[] | null {
  const taxon = genome?.taxon_id ?? undefined;
  const species = genome?.species_name ?? undefined;
  const params =
    taxon !== undefined
      ? { taxon_id: taxon, pageSize: 1000 }
      : {
          q: species,
          searchFields: ['species'],
          operator: 'eq' as const,
          pageSize: 1000,
        };
  const enabled = !!genome && (taxon !== undefined || !!species);
  const list = useGenomes(params, { enabled });

  return useMemo(() => {
    if (!genome) return null;
    if (!enabled) return [];
    if (!list.data) return list.error ? [] : null;
    return rankCompareCandidates(genome, list.data.results);
  }, [genome, enabled, list.data, list.error]);
}
