/**
 * `/tree` — the radial tree of life, with all of its chrome as real HTML.
 *
 * The tree's controls live here as HTML, not as an invisible SVG overlay
 * inside `Tree.tsx`: the filter, the level select, zoom, the interaction mode, full screen,
 * and the species detail panel. Nothing is revealed on hover, and full screen
 * is never a side effect of typing.
 */

import { useCallback, useDeferredValue, useEffect, useRef, useState } from 'react';
import { Link, useSearchParams } from 'react-router-dom';
import { MiniHero } from '../components/layout/MiniHero';
import { Tree } from '../components/tree/Tree';
import type { TreeHandle } from '../components/tree/Tree';
import { resolveSpecies } from '../components/tree/taxa';
import { useSpeciesIndex } from '../hooks/queries/useSpecies';
import { TAXONOMIC_LEVELS } from '../types/taxonomy';
import type { TaxonomicLevel } from '../types/taxonomy';
import { cn } from '../utils/cn';
import { normalizeSpecies, speciesGenomesHref } from '../utils/species';

const INTERACTION_MODES = ['inspect', 'pan'] as const;

/** The selected tip. In the URL so a tree view is a link you can send. */
const SPECIES_PARAM = 'species';

export function TreeExplorerPage() {
  const [filter, setFilter] = useState('');
  const deferredFilter = useDeferredValue(filter);
  const [level, setLevel] = useState<TaxonomicLevel>('class');
  const [interaction, setInteraction] = useState<'inspect' | 'pan'>('inspect');
  const [fullscreen, setFullscreen] = useState(false);
  const [matches, setMatches] = useState({ count: 0, suppressed: false });

  const [params, setParams] = useSearchParams();
  const requested = params.get(SPECIES_PARAM);
  // Fold to the taxonomy's own spelling: an inbound link from /species carries
  // the API's "Homo sapiens", and the tree's tip is named "Homo Sapiens".
  const selected = requested ? resolveSpecies(requested) : null;

  const setSelected = useCallback(
    (species: string | null) => {
      setParams(
        (previous) => {
          const next = new URLSearchParams(previous);
          if (species) next.set(SPECIES_PARAM, species);
          else next.delete(SPECIES_PARAM);
          return next;
        },
        // Replace, not push: clicking through six tips must not bury the page
        // the visitor arrived from under six history entries.
        { replace: true },
      );
    },
    [setParams],
  );

  const treeRef = useRef<TreeHandle>(null);
  const species = useSpeciesIndex();
  // Keyed on the normalized name, which is what fixes the "Homo Sapiens" vs
  // "Homo sapiens" miss that reported no genomes for humans. `taxa.json`
  // capitalizes the epithet and the API does not, so an exact lookup misses on
  // every such tip.
  const info = selected ? species.data?.get(normalizeSpecies(selected)) : undefined;

  // Same object back when nothing changed, so a redraw costs no render.
  const handleMatchStats = useCallback((count: number, suppressed: boolean) => {
    setMatches((prev) =>
      prev.count === count && prev.suppressed === suppressed ? prev : { count, suppressed },
    );
  }, []);

  useEffect(() => {
    if (!fullscreen) return;
    const onKey = (event: KeyboardEvent) => {
      if (event.key === 'Escape') setFullscreen(false);
    };
    window.addEventListener('keydown', onKey);
    return () => window.removeEventListener('keydown', onKey);
  }, [fullscreen]);

  return (
    <div className="flex flex-col gap-4">
      {/* Rendered in full screen too, where `.rg-tree__panel--fullscreen` covers
          it completely: hiding it would also unmount the hook that owns the
          document title. */}
      <MiniHero
        title="Tree of life"
        lede={
          <>
            Every species refgenie knows about, arranged by taxonomy. Hover a tip to read its
            name, click to select it, then jump to its genomes or to its row in the species
            table.
          </>
        }
      />

      {/* Always visible. Never inside the SVG, never revealed on hover. */}
      <div
        className={cn(
          'rg-tree__toolbar flex flex-wrap items-center gap-2',
          fullscreen && 'rg-tree__toolbar--fullscreen',
        )}
      >
        <label className="sr-only" htmlFor="tree-filter">
          Filter species
        </label>
        <input
          id="tree-filter"
          className="rg-field__input rg-tree__filter"
          type="search"
          placeholder="Filter species…"
          value={filter}
          onChange={(event) => setFilter(event.target.value)}
        />

        <label className="sr-only" htmlFor="tree-level">
          Group by taxonomic level
        </label>
        <select
          id="tree-level"
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

        <div className="flex gap-1" role="group" aria-label="Interaction mode">
          {INTERACTION_MODES.map((mode) => (
            <button
              key={mode}
              type="button"
              className={cn('rg-btn rg-btn--sm', interaction === mode && 'rg-btn--primary')}
              aria-pressed={interaction === mode}
              onClick={() => setInteraction(mode)}
            >
              {mode === 'inspect' ? 'Inspect' : 'Pan'}
            </button>
          ))}
        </div>

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

        <button
          type="button"
          className="rg-btn rg-btn--sm rg-tree__toolbar-end"
          onClick={() => setFullscreen((value) => !value)}
        >
          {fullscreen ? 'Exit full screen' : 'Full screen'}
        </button>
      </div>

      {filter.trim() !== '' && (
        <p className="rg-muted text-sm" role="status">
          {matches.count} matching species
          {matches.suppressed && ' — narrow the filter to see their names'}
        </p>
      )}

      <div className={cn('rg-tree__panel', fullscreen && 'rg-tree__panel--fullscreen')}>
        <Tree
          ref={treeRef}
          filter={deferredFilter}
          level={level}
          interaction={interaction}
          selected={selected}
          onSelect={setSelected}
          onMatchStats={handleMatchStats}
        />

        {selected && (
          <aside className="rg-tree__detail rg-card">
            <div className="rg-card__body flex flex-col gap-2">
              <div className="flex items-start justify-between gap-2">
                <h2 className="font-semibold">{selected}</h2>
                <button
                  type="button"
                  className="rg-btn rg-btn--bare rg-btn--sm"
                  onClick={() => setSelected(null)}
                >
                  Clear
                </button>
              </div>
              {info ? (
                <>
                  {info.commonNames.length > 0 && (
                    <p className="rg-muted text-sm">{info.commonNames.join(', ')}</p>
                  )}
                  <p className="rg-muted text-sm">
                    {info.genomes} genome{info.genomes === 1 ? '' : 's'} · {info.assets} asset
                    {info.assets === 1 ? '' : 's'}
                  </p>
                  <Link
                    className="rg-btn rg-btn--primary rg-btn--sm"
                    to={speciesGenomesHref(info) ?? '/genomes'}
                  >
                    View genomes
                  </Link>
                  <Link
                    className="rg-btn rg-btn--sm"
                    to={`/species?q=${encodeURIComponent(info.speciesName)}`}
                  >
                    Species table
                  </Link>
                </>
              ) : (
                <p className="rg-muted text-sm">No genomes for this species on this instance.</p>
              )}
            </div>
          </aside>
        )}
      </div>
    </div>
  );
}
