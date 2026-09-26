/**
 * Genome picker for the build form.
 *
 * Submits the ALIAS, because `rgc.build.run` resolves the alias itself;
 * the digest travels alongside only so the job's target key is stable.
 *
 * A plain `<select>` over the whole local list, matching the recipe picker, not
 * a text input with a `datalist`: that is a real combobox to the platform, but
 * on screen it is a text box whose suggestions look like tooltips, so nobody
 * reads it as "pick one of these". Search only pays off for a catalog too
 * large to enumerate; a dash instance holds tens of genomes, and the
 * server-side `q` + `searchFields: ['aliases']` filter is there for whoever
 * needs it.
 */

import { useGenomeIndexResource } from '../../hooks/queries/useGenomes';
import { FormField } from '../common/FormField';
import type { GenomeResponse } from '../../types/api';

export interface GenomeSelectProps {
  value: string;
  onChange: (alias: string, digest: string | null) => void;
  error?: string;
  /** Rendered as the call to action when there are no genomes at all. */
  onInitGenome?: () => void;
}

interface GenomeOption {
  alias: string;
  digest: string;
}

/**
 * Module-level, because `useGenomeIndexResource` memoizes on the selector's
 * identity — a closure defined in the component body would rebuild the array
 * on every render.
 */
function toOptions(genomes: readonly GenomeResponse[]): GenomeOption[] {
  return genomes
    .flatMap((genome) =>
      (genome.aliases.length ? genome.aliases : [genome.digest]).map((alias) => ({
        alias,
        digest: genome.digest,
      })),
    )
    .sort((a, b) => a.alias.localeCompare(b.alias));
}

export function GenomeSelect({ value, onChange, error, onInitGenome }: GenomeSelectProps) {
  const genomes = useGenomeIndexResource(toOptions);
  const options = genomes.data ?? [];

  if (!genomes.isPending && options.length === 0) {
    return (
      <div className="rg-field">
        <div className="rg-field__label-cell">
          <p className="rg-field__label">Genome</p>
        </div>
        <div className="rg-field__control flex flex-col items-start gap-2">
          <p className="rg-muted text-sm">
            No genomes yet. A build needs something to build against.
          </p>
          {onInitGenome && (
            <button type="button" className="rg-btn rg-btn--primary" onClick={onInitGenome}>
              Initialize a genome
            </button>
          )}
        </div>
      </div>
    );
  }

  // `?genome=` can name something this instance does not have. Dropping it
  // would show an empty control while the form still holds the value.
  const unknown = value !== '' && !options.some((option) => option.alias === value);

  return (
    <FormField
      htmlFor="build-genome"
      label="Genome"
      required
      error={error}
      hint="The build resolves the alias to a digest server-side."
    >
      <select
        id="build-genome"
        className="rg-field__input"
        value={value}
        onChange={(event) => {
          const alias = event.target.value;
          const match = options.find((option) => option.alias === alias);
          onChange(alias, match?.digest ?? null);
        }}
      >
        <option value="">{genomes.isPending ? 'Loading genomes…' : 'Choose a genome…'}</option>
        {unknown && <option value={value}>{value} (not on this instance)</option>}
        {options.map((option) => (
          <option key={`${option.digest}-${option.alias}`} value={option.alias}>
            {option.alias}
          </option>
        ))}
      </select>
    </FormField>
  );
}
