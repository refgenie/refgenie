import { Copyable } from '../common/Copyable';
import { ExternalLink } from '../common/ExternalLink';
import { accessionUrl } from '../../utils/links';
import type { AliasEntry } from '../../types/api';

export interface GenomeAliasesProps {
  /** The genome's refgenie aliases. */
  aliases: string[];
  /** The refget store's aliases for this genome's sequence collection. */
  collectionAliases?: AliasEntry[];
}

interface AliasRow {
  alias: string;
  sources: string[];
}

/** One row per alias name, listing every place it comes from. */
function aliasRows(
  aliases: string[],
  collectionAliases: AliasEntry[] = [],
): AliasRow[] {
  const rows = new Map<string, AliasRow>();
  const add = (alias: string, source: string) => {
    const row = rows.get(alias) ?? { alias, sources: [] };
    if (!row.sources.includes(source)) row.sources.push(source);
    rows.set(alias, row);
  };
  for (const alias of aliases) add(alias, 'refgenie');
  for (const entry of collectionAliases) {
    add(entry.alias, `refget store, ${entry.namespace}`);
  }
  return [...rows.values()];
}

/** The Details tab's alias list: one alias a row with its copy button, and where it comes from. */
export function GenomeAliases({ aliases, collectionAliases }: GenomeAliasesProps) {
  const rows = aliasRows(aliases, collectionAliases);
  return (
    <section id="aliases" className="rg-flat-section">
      <h2>Aliases</h2>
      {rows.length === 0 ? (
        <p className="rg-muted">This genome has no aliases.</p>
      ) : (
        <div className="rg-table__wrap">
          <table className="rg-table rg-table--compact rg-alias-table">
            <caption className="sr-only">Aliases</caption>
            <thead className="rg-table__head">
              <tr>
                <th scope="col" className="rg-table__cell">
                  Alias
                </th>
                <th scope="col" className="rg-table__cell">
                  From
                </th>
              </tr>
            </thead>
            <tbody>
              {rows.map((row) => (
                <tr key={row.alias} className="rg-table__row">
                  <td className="rg-table__cell rg-alias-table__alias">
                    <AliasValue alias={row.alias} />
                  </td>
                  <td className="rg-table__cell rg-muted">{row.sources.join('; ')}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}
    </section>
  );
}

/** An alias with its copy button; an INSDC or RefSeq accession also links to NCBI. */
function AliasValue({ alias }: { alias: string }) {
  const href = accessionUrl(alias);
  return (
    <Copyable value={alias} label={`Copy ${alias}`} mono>
      {href ? (
        <ExternalLink href={href} title="View this assembly at NCBI">
          {alias}
        </ExternalLink>
      ) : (
        alias
      )}
    </Copyable>
  );
}
