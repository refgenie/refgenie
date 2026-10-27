/** The shape of one entry in a recipe's `input_params` / `input_files` /
 *  `input_assets` blob. Untyped on the server (`InputEntities`), so every
 *  field is optional and every read is guarded. */
export interface EntrySpec {
  description?: string;
  default?: unknown;
  asset_class?: string;
  [key: string]: unknown;
}

export function specEntries(
  source: Record<string, Record<string, unknown>> | null | undefined,
): Array<[string, EntrySpec]> {
  return Object.entries(source ?? {}) as Array<[string, EntrySpec]>;
}
