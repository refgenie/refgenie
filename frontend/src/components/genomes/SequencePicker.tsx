import { useId, useMemo, useState } from 'react';
import { primaryChromosomes } from '../../utils/fastaDownload';

export interface SequencePickerProps {
  /** All sequence names with lengths, in file order. */
  sequences: Array<{ name: string; length: number }>;
  selected: Set<string>;
  onChange: (selected: Set<string>) => void;
  /** Names that cannot be chosen, each with a short reason shown beside it. */
  unavailable?: ReadonlyMap<string, string>;
}

const NONE: ReadonlyMap<string, string> = new Map();

/** Quick presets and a filterable checklist for choosing sequences to export. */
export function SequencePicker({
  sequences,
  selected,
  onChange,
  unavailable = NONE,
}: SequencePickerProps) {
  const [filter, setFilter] = useState('');
  const filterId = useId();
  const names = useMemo(
    () => sequences.map((s) => s.name).filter((name) => !unavailable.has(name)),
    [sequences, unavailable],
  );
  const primary = useMemo(() => primaryChromosomes(names), [names]);
  const needle = filter.trim().toLowerCase();
  const shown = needle
    ? sequences.filter((s) => s.name.toLowerCase().includes(needle))
    : sequences;

  const toggle = (name: string) => {
    const next = new Set(selected);
    if (next.has(name)) next.delete(name);
    else next.add(name);
    onChange(next);
  };

  const setShown = (on: boolean) => {
    const next = new Set(selected);
    for (const s of shown) {
      if (on) {
        if (!unavailable.has(s.name)) next.add(s.name);
      } else next.delete(s.name);
    }
    onChange(next);
  };

  return (
    <div className="rg-seq-picker">
      <div className="flex flex-wrap items-center gap-2">
        {primary && (
          <button
            type="button"
            className="rg-btn rg-btn--sm"
            onClick={() => onChange(new Set(primary))}
          >
            Primary chromosomes ({primary.length})
          </button>
        )}
        <button
          type="button"
          className="rg-btn rg-btn--sm"
          onClick={() => setShown(true)}
        >
          {needle ? 'Select shown' : 'Select all'}
        </button>
        <button
          type="button"
          className="rg-btn rg-btn--sm"
          onClick={() => setShown(false)}
        >
          {needle ? 'Clear shown' : 'Clear'}
        </button>
        <label htmlFor={filterId} className="sr-only">
          Filter sequences
        </label>
        <input
          id={filterId}
          type="search"
          className="rg-field__input rg-seq-picker__filter"
          placeholder="Filter by name…"
          value={filter}
          onChange={(event) => setFilter(event.target.value)}
        />
      </div>
      <ul className="rg-seq-picker__list" aria-label="Sequences to download">
        {shown.map((s) => {
          const reason = unavailable.get(s.name);
          return (
            <li key={s.name}>
              <label className="rg-seq-picker__item">
                <input
                  type="checkbox"
                  checked={!reason && selected.has(s.name)}
                  disabled={!!reason}
                  onChange={() => toggle(s.name)}
                />
                <span>{s.name}</span>
                <span className="rg-muted">
                  {s.length.toLocaleString()} bp{reason && ` · ${reason}`}
                </span>
              </label>
            </li>
          );
        })}
        {shown.length === 0 && <li className="rg-muted">No sequence names match.</li>}
      </ul>
    </div>
  );
}
