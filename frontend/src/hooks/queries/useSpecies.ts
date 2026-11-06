import { useGenomeIndexResource } from './useGenomes';
import { buildSpeciesIndex } from '../../utils/species';

/**
 * Species, grouped from the cached genome list. Not from `/v4/species/summary`,
 * which lacks common names, drops species with no assets, and is absent on a local dash.
 */
export function useSpeciesIndex() {
  return useGenomeIndexResource(buildSpeciesIndex);
}
