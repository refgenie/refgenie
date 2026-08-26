/**
 * Genome picker for the build form.
 *
 * Submits the ALIAS, because `facade.build_asset` resolves the alias itself;
 * the digest travels alongside only so the job's target key is stable. A
 * `datalist` rather than a custom listbox: it is a real combobox to the
 * platform, keyboard and screen reader included, with no widget to maintain.
 */

import { useEffect, useMemo, useState } from 'react';
import { useGenomes } from '../../hooks/queries/useGenomes';
import { FormField } from '../common/FormField';
import { formatDigest } from '../../utils/format';

export interface GenomeSelectProps {
  value: string;
  onChange: (alias: string, digest: string | null) => void;
  error?: string;
  /** Rendered as the call to action when there are no genomes at all. */
  onInitGenome?: () => void;
}

const SEARCH_LIMIT = 25;

export function GenomeSelect({ value, onChange, error, onInitGenome }: GenomeSelectProps) {
  const [typed, setTyped] = useState(value);
  const [query, setQuery] = useState('');

  useEffect(() => setTyped(value), [value]);

  // Debounced so a fast typist does not issue a request per keystroke.
  useEffect(() => {
    const timer = setTimeout(() => setQuery(typed.trim()), 250);
    return () => clearTimeout(timer);
  }, [typed]);

  const genomes = useGenomes({
    q: query || undefined,
    // 'aliases' is on the server's allowlist; anything else is a 422.
    searchFields: query ? ['aliases'] : undefined,
    limit: SEARCH_LIMIT,
  });

  const options = useMemo(
    () =>
      (genomes.data?.items ?? []).flatMap((genome) =>
        (genome.aliases.length ? genome.aliases : [genome.digest]).map((alias) => ({
          alias,
          digest: genome.digest,
        })),
      ),
    [genomes.data],
  );

  const empty = !genomes.isPending && !query && options.length === 0;

  if (empty) {
    return (
      <div className="rg-field">
        <p className="rg-field__label">Genome</p>
        <p className="rg-muted text-sm mb-2">
          No genomes yet. A build needs something to build against.
        </p>
        {onInitGenome && (
          <button type="button" className="rg-btn rg-btn--primary" onClick={onInitGenome}>
            Initialize a genome
          </button>
        )}
      </div>
    );
  }

  return (
    <FormField
      htmlFor="build-genome"
      label="Genome"
      required
      error={error}
      hint="Type an alias. The build resolves it to a digest server-side."
    >
      <input
        id="build-genome"
        className="rg-field__input"
        type="text"
        list="build-genome-options"
        value={typed}
        autoComplete="off"
        placeholder="hg38"
        onChange={(event) => {
          const next = event.target.value;
          setTyped(next);
          const match = options.find((option) => option.alias === next);
          onChange(next, match?.digest ?? null);
        }}
      />
      <datalist id="build-genome-options">
        {options.map((option) => (
          <option key={`${option.digest}-${option.alias}`} value={option.alias}>
            {formatDigest(option.digest)}
          </option>
        ))}
      </datalist>
    </FormField>
  );
}
