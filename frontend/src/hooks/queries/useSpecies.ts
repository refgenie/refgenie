import { useGenomeIndexResource } from './useGenomes';
import { buildSpeciesIndex } from '../../utils/species';

/**
 * Species, grouped from the cached genome list.
 *
 * Same cache entry as `useGenomeIndex()`, different `select` — see
 * `useGenomeIndexResource`. There is no `/v4/species` endpoint and this hook
 * does not call `/v4/species/summary`: that endpoint carries no common name,
 * inner-joins away every species whose genomes have no assets, and is not
 * mounted on a local dash.
 */
export function useSpeciesIndex() {
  return useGenomeIndexResource(buildSpeciesIndex);
}
