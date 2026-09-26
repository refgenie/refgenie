/**
 * Slot marker for command affordances.
 *
 * The slots that DO have controls now render them directly (pull on remote
 * asset rows, delete on genome/asset headers and rows, set-default in an asset
 * group). The three that remain empty are empty on purpose:
 *  - `asset-class` and `recipe` — registration takes a server-side path or URL
 *    from an HTTP body, which is an arbitrary local-file read plus SSRF, so
 *    those writes are out of v1 and remain CLI operations;
 *  - `remote-genome` — there is no whole-genome pull verb.
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
