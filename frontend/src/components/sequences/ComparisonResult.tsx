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
}

/** A seqcol comparison: a verdict badge, plain sentences, and the per-attribute table. */
export function ComparisonResult({
  comparison,
  caption,
  otherName = 'the other genome',
  otherColumn = 'In other',
}: ComparisonResultProps) {
  const verdict = interpretCoordinateSystem(comparison, otherName);
  const counts = comparison.array_elements;
  return (
    <>
      <div className="flex flex-col gap-2">
        <p className="flex items-center gap-2 flex-wrap">
          <Badge variant={LEVEL_BADGE[verdict.level].variant}>
            {LEVEL_BADGE[verdict.level].label}
          </Badge>
          <span className="text-sm">{caption}</span>
        </p>
        {verdict.messages.map((message) => (
          <p key={message}>{message}</p>
        ))}
      </div>
      <div className="rg-table__wrap">
        <table className="rg-table rg-table--compact">
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
                <td className="rg-table__cell">{attribute}</td>
                <td className="rg-table__cell rg-table__cell--numeric">
                  {counts.a_count[attribute] ?? 'NA'}
                </td>
                <td className="rg-table__cell rg-table__cell--numeric">
                  {counts.b_count[attribute] ?? 'NA'}
                </td>
                <td className="rg-table__cell rg-table__cell--numeric">
                  {counts.a_and_b_count[attribute] ?? 'NA'}
                </td>
                <td className="rg-table__cell">
                  {orderLabel(counts.a_and_b_same_order[attribute])}
                </td>
                <td className="rg-table__cell">
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
