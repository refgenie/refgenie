// Ported from repos/refget/frontend/src/pages/ComparisonView.tsx @ c0d3e09; keep the two in step.

import { useState } from 'react';
import type { FormEvent } from 'react';
import { useGenomeComparison } from '../../hooks/queries/useGenomeSequences';
import { FormField } from '../common/FormField';
import { ErrorState, LoadingState } from '../common/states';
import { ScopedLink } from '../common/ScopedLink';
import { ComparisonResult } from './ComparisonResult';
import { findCandidate } from '../../utils/compareCandidates';
import type { CompareCandidate } from '../../utils/compareCandidates';

/** Pick another genome and read how its sequence collection relates to this one. */
export function CompareGenomes({
  digest,
  candidates,
}: {
  digest: string;
  candidates: CompareCandidate[];
}) {
  const [selected, setSelected] = useState('');
  const [text, setText] = useState('');
  const [error, setError] = useState<string>();
  const comparison = useGenomeComparison(digest, selected || undefined);
  const options = candidates.filter((c) => c.digest !== digest);
  const other = options.find((c) => c.digest === selected);

  // Enter or the button: the genome named, else the first one the text matches.
  const submit = (event: FormEvent) => {
    event.preventDefault();
    const match = findCandidate(options, text);
    if (match) {
      setSelected(match.digest);
      setText(match.label);
      setError(undefined);
    } else {
      setError(
        text.trim()
          ? `No genome of this species matches "${text.trim()}". Choose one from the list.`
          : 'No genome selected. Choose one from the list.',
      );
    }
  };

  return (
    <section className="rg-flat-section">
      <h2>Compare with another genome of this species</h2>
      <div className="flex flex-col gap-4">
        {options.length === 0 ? (
          <p className="rg-muted text-sm">
            No other genomes of this species to compare with.
          </p>
        ) : (
          <form onSubmit={submit}>
            <FormField htmlFor="compare-genome" label="Other genome" error={error}>
              <div className="flex gap-2">
                <input
                  id="compare-genome"
                  className="rg-field__input"
                  list="compare-genome-list"
                  autoComplete="off"
                  placeholder={`Type to filter ${options.length} genomes…`}
                  value={text}
                  // Selected on focus, so typing replaces the current choice.
                  onFocus={(event) => event.target.select()}
                  onChange={(event) => {
                    setText(event.target.value);
                    setError(undefined);
                    const match = options.find((c) => c.label === event.target.value);
                    if (match) setSelected(match.digest);
                  }}
                />
                <button type="submit" className="rg-btn">
                  Compare
                </button>
              </div>
              <datalist id="compare-genome-list">
                {options.map((c) => (
                  <option key={c.digest} value={c.label} />
                ))}
              </datalist>
            </FormField>
          </form>
        )}

        {selected && comparison.isPending && <LoadingState label="Comparing" />}
        {selected && !!comparison.error && <ErrorState error={comparison.error} />}

        {selected && comparison.data && (
          <ComparisonResult
            comparison={comparison.data}
            caption={
              <>
                Compared with{' '}
                <ScopedLink
                  className="rg-link"
                  to={`/genomes/${encodeURIComponent(selected)}`}
                >
                  {other?.label ?? selected}
                </ScopedLink>
              </>
            }
          />
        )}
      </div>
    </section>
  );
}
