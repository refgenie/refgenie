/**
 * The taxonomy tree cut down to a collection's species. Default export so the
 * collection page can `React.lazy` it and keep d3 and `taxa.json` out of the
 * main bundle.
 */

import { useCallback, useMemo, useRef, useState } from 'react';
import { Link } from 'react-router-dom';
import { ScopedLink } from '../common/ScopedLink';
import { Tree } from '../tree/Tree';
import type { TreeHandle } from '../tree/Tree';
import { buildTreeForSpecies } from '../tree/taxa';
import { TAXONOMIC_LEVELS } from '../../types/taxonomy';
import type { TaxonomicLevel } from '../../types/taxonomy';
import { preferredAlias } from '../../utils/genomes';
import { formatDigest } from '../../utils/format';
import { normalizeSpecies } from '../../utils/species';
import type { GenomeResponse } from '../../types/api';

const noop = () => {};

export default function CollectionTree({ genomes }: { genomes: GenomeResponse[] }) {
  const [level, setLevel] = useState<TaxonomicLevel>('class');
  const [selected, setSelected] = useState<string | null>(null);
  const treeRef = useRef<TreeHandle>(null);

  const { tree, missing } = useMemo(
    () =>
      buildTreeForSpecies(
        genomes.map((g) => g.species_name).filter((name): name is string => !!name),
      ),
    [genomes],
  );

  const selectedGenomes = useMemo(() => {
    if (!selected) return [];
    const key = normalizeSpecies(selected);
    return genomes.filter((g) => normalizeSpecies(g.species_name) === key);
  }, [genomes, selected]);

  const onSelect = useCallback((species: string | null) => setSelected(species), []);

  return (
    <div className="flex flex-col gap-2">
      <div className="rg-tree__toolbar flex flex-wrap items-center gap-2">
        <label className="sr-only" htmlFor="collection-tree-level">
          Group by taxonomic level
        </label>
        <select
          id="collection-tree-level"
          className="rg-field__input w-auto"
          value={level}
          onChange={(event) => setLevel(event.target.value as TaxonomicLevel)}
        >
          {TAXONOMIC_LEVELS.map((option) => (
            <option key={option} value={option}>
              Group by {option}
            </option>
          ))}
        </select>
        <div className="flex gap-1" role="group" aria-label="Zoom">
          <button
            type="button"
            className="rg-btn rg-btn--sm"
            onClick={() => treeRef.current?.zoomIn()}
          >
            Zoom in
          </button>
          <button
            type="button"
            className="rg-btn rg-btn--sm"
            onClick={() => treeRef.current?.zoomOut()}
          >
            Zoom out
          </button>
          <button
            type="button"
            className="rg-btn rg-btn--sm"
            onClick={() => treeRef.current?.reset()}
          >
            Reset
          </button>
        </div>
      </div>

      {missing.length > 0 && (
        <p className="rg-muted text-sm">
          {missing.length} species in this collection{' '}
          {missing.length === 1 ? 'is' : 'are'} not in the bundled taxonomy and{' '}
          {missing.length === 1 ? 'is' : 'are'} not drawn.
        </p>
      )}

      <div className="rg-tree__panel">
        <Tree
          ref={treeRef}
          data={tree}
          filter=""
          level={level}
          interaction="inspect"
          selected={selected}
          onSelect={onSelect}
          onMatchStats={noop}
        />
        {selected && (
          <aside className="rg-tree__detail">
            <div className="flex flex-col gap-2">
              <div className="flex items-start justify-between gap-2">
                <h3 className="font-semibold">{selected}</h3>
                <button
                  type="button"
                  className="rg-btn rg-btn--bare rg-btn--sm"
                  onClick={() => setSelected(null)}
                >
                  Clear
                </button>
              </div>
              <p className="rg-muted text-sm">
                {selectedGenomes.length} genome{selectedGenomes.length === 1 ? '' : 's'}{' '}
                in this collection
              </p>
              <ul className="flex flex-col gap-1 text-sm">
                {selectedGenomes.map((g) => (
                  <li key={g.digest}>
                    <ScopedLink className="rg-link" to={`/genomes/${g.digest}`}>
                      {preferredAlias(g.aliases) ?? formatDigest(g.digest)}
                    </ScopedLink>
                  </li>
                ))}
              </ul>
              <Link
                className="rg-btn rg-btn--sm"
                to={`/tree?species=${encodeURIComponent(selected)}`}
              >
                Open in the full tree
              </Link>
            </div>
          </aside>
        )}
      </div>
    </div>
  );
}
