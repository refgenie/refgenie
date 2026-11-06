/**
 * Slot marker for command controls. The `asset-class`, `recipe` and
 * `remote-genome` slots stay empty: registering classes and recipes over HTTP
 * would allow local file reads and SSRF, and there is no whole-genome pull.
 */

export type ActionSlot =
  | 'genome'
  | 'genome-asset'
  | 'asset'
  | 'asset-class'
  | 'recipe'
  | 'remote-genome'
  | 'remote-asset';

export interface ActionBarProps {
  slot: ActionSlot;
  context?: Record<string, string | number | null | undefined>;
}

export function ActionBar(_props: ActionBarProps) {
  return null;
}
