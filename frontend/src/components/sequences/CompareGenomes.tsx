// Ported from repos/refget/frontend/src/pages/ComparisonView.tsx @ c0d3e09; keep the two in step.

import { FormField } from '../common/FormField';
import { ErrorState, LoadingState } from '../common/states';
import { ScopedLink } from '../common/ScopedLink';
import { ComparisonResult } from './ComparisonResult';
import type { SeqcolComparison } from './types';
import type { CompareCandidate } from '../../utils/compareCandidates';

export interface CompareGenomesProps {
  digest: string;
  /** This genome's name, for the summary sentence; "this genome" if none. */
  name?: string;
  candidates: CompareCandidate[];
  comparison?: SeqcolComparison;
  loading: boolean;
  error: unknown;
  onSelect: (other: string) => void;
  selected?: string;
}

/** Pick another genome and read how its sequence collection relates to this one. */
export function CompareGenomes({
  digest,
  name,
  candidates,
  comparison,
  loading,
  error,
  onSelect,
  selected,
}: CompareGenomesProps) {
  const options = candidates.filter((c) => c.digest !== digest);
  const other = options.find((c) => c.digest === selected);

  return (
    <section className="rg-flat-section">
      <h2>Compare with another genome of this species</h2>
      <div className="flex flex-col gap-4">
        {options.length === 0 ? (
          <p className="rg-muted text-sm">
            No other genomes of this species to compare with.
          </p>
        ) : (
          <FormField htmlFor="compare-genome" label="Other genome">
            <input
              id="compare-genome"
              className="rg-field__input"
              list="compare-genome-list"
              autoComplete="off"
              placeholder={`Type to filter ${options.length} genomes…`}
              defaultValue={other?.label}
              // Selected on focus, so typing replaces the current choice.
              onFocus={(event) => event.target.select()}
              onChange={(event) => {
                const match = options.find((c) => c.label === event.target.value);
                if (match) onSelect(match.digest);
              }}
            />
            <datalist id="compare-genome-list">
              {options.map((c) => (
                <option key={c.digest} value={c.label} />
              ))}
            </datalist>
          </FormField>
        )}

        {selected && loading && <LoadingState label="Comparing" />}
        {selected && !!error && <ErrorState error={error} />}

        {selected && comparison && (
          <ComparisonResult
            comparison={comparison}
            names={{ a: name ?? 'this genome', b: other?.name ?? 'the other genome' }}
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
