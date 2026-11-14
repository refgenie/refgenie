// Ported from repos/refget/frontend/src/pages/ComparisonView.tsx @ c0d3e09; keep the two in step.

import type { ReactNode } from 'react';
import { Badge } from '../common/Badge';
import {
  interpretAttribute,
  interpretCoordinateSystem,
} from '../../utils/seqcolCompare';
import type { BadgeVariant } from '../common/Badge';
import type { CompareLevel } from '../../utils/seqcolCompare';
import type { SeqcolComparison } from './types';

const LEVEL_BADGE: Record<CompareLevel, { variant: BadgeVariant; label: string }> = {
  same: { variant: 'archive', label: 'Same coordinates' },
  subset: { variant: 'file', label: 'Subset' },
  overlap: { variant: 'default', label: 'Partial overlap' },
  disjoint: { variant: 'default', label: 'Disjoint' },
  swap: { variant: 'default', label: 'Name swap' },
};

function orderLabel(value: boolean | null | undefined): string {
  if (value === true) return 'Yes';
  if (value === false) return 'No';
  return 'NA';
}

export interface ComparisonResultProps {
  comparison: SeqcolComparison;
  /** Beside the badge, e.g. "Compared with <link>". */
  caption: ReactNode;
  /** The other side in sentences, e.g. "the other genome", "your file". */
  otherName?: string;
  /** The other side's column heading. */
  otherColumn?: string;
  /** Both genomes' names, for the server's summary and the rename table. */
  names?: { a: string; b: string };
}

/** The server's summary names each genome by its digest; swap in the names shown here. */
function nameSummary(
  summary: string,
  digests: { a: string; b: string } | undefined,
  names: { a: string; b: string },
): string {
  let text = summary;
  if (digests?.a) text = text.split(digests.a).join(names.a);
  if (digests?.b) text = text.split(digests.b).join(names.b);
  return text;
}

/** A's name -> B's name for each sequence B has under another name. */
function Renames({
  map,
  names,
}: {
  map: Record<string, string>;
  names: { a: string; b: string };
}) {
  const pairs = Object.entries(map);
  if (pairs.length === 0) return null;
  return (
    <div className="rg-table__wrap rg-compare__renames">
      <table className="rg-table rg-table--compact">
        <caption className="rg-compare__renames-caption">
          {pairs.length === 1 ? '1 sequence has' : `${pairs.length} sequences have`} another
          name
        </caption>
        <thead className="rg-table__head">
          <tr>
            <th scope="col" className="rg-table__cell">
              Name in {names.a}
            </th>
            <th scope="col" className="rg-table__cell">
              Name in {names.b}
            </th>
          </tr>
        </thead>
        <tbody>
          {pairs.map(([a, b]) => (
            <tr key={a} className="rg-table__row">
              <td className="rg-table__cell">
                <code>{a}</code>
              </td>
              <td className="rg-table__cell">
                <code>{b}</code>
              </td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

/** A seqcol comparison: a verdict badge, plain sentences, and the per-attribute table. */
export function ComparisonResult({
  comparison,
  caption,
  otherName = 'the other genome',
  otherColumn = 'In other',
  names = { a: 'this genome', b: otherName },
}: ComparisonResultProps) {
  const verdict = interpretCoordinateSystem(comparison, otherName);
  const counts = comparison.array_elements;
  const summary = comparison.summary
    ? nameSummary(comparison.summary, comparison.digests, names)
    : null;
  return (
    <>
      <div className="flex flex-col gap-2">
        <p className="flex items-center gap-2 flex-wrap">
          <Badge variant={LEVEL_BADGE[verdict.level].variant}>
            {LEVEL_BADGE[verdict.level].label}
          </Badge>
          <span className="text-sm">{caption}</span>
        </p>
        {summary && <p className="rg-compare__summary">{summary}</p>}
        {comparison.name_map && <Renames map={comparison.name_map} names={names} />}
        {verdict.messages.map((message) => (
          <p key={message} className={summary ? 'rg-muted text-sm' : undefined}>
            {message}
          </p>
        ))}
      </div>
      <div className="rg-table__wrap">
        <table className="rg-table rg-table--compact rg-table--stack">
          <caption className="sr-only">Comparison by attribute</caption>
          <thead className="rg-table__head">
            <tr>
              <th scope="col" className="rg-table__cell">
                Attribute
              </th>
              <th scope="col" className="rg-table__cell rg-table__cell--numeric">
                In this genome
              </th>
              <th scope="col" className="rg-table__cell rg-table__cell--numeric">
                {otherColumn}
              </th>
              <th scope="col" className="rg-table__cell rg-table__cell--numeric">
                Shared
              </th>
              <th scope="col" className="rg-table__cell">
                Same order
              </th>
              <th scope="col" className="rg-table__cell">
                Interpretation
              </th>
            </tr>
          </thead>
          <tbody>
            {comparison.attributes.a_and_b.map((attribute) => (
              <tr key={attribute} className="rg-table__row">
                <td className="rg-table__cell rg-table__cell--wrap">{attribute}</td>
                <td className="rg-table__cell rg-table__cell--numeric" data-label="In this genome">
                  {counts.a_count[attribute] ?? 'NA'}
                </td>
                <td className="rg-table__cell rg-table__cell--numeric" data-label={otherColumn}>
                  {counts.b_count[attribute] ?? 'NA'}
                </td>
                <td className="rg-table__cell rg-table__cell--numeric" data-label="Shared">
                  {counts.a_and_b_count[attribute] ?? 'NA'}
                </td>
                <td className="rg-table__cell" data-label="Same order">
                  {orderLabel(counts.a_and_b_same_order[attribute])}
                </td>
                <td className="rg-table__cell" data-label="Interpretation">
                  {interpretAttribute(comparison, attribute, otherName)}
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </>
  );
}
