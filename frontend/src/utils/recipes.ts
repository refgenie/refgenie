/** One entry of a recipe's `input_params`, `input_files` or `input_assets`; untyped on the server. */
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
