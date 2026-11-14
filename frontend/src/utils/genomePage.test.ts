import { describe, expect, it } from 'vitest';
import { fastaPullSize } from './genomePage';
import type { AssetResponse } from '../types/api';

const asset = (over: Partial<AssetResponse>): AssetResponse => ({
  digest: 'd',
  name: 'default',
  description: null,
  recipe_id: null,
  asset_group_id: null,
  size: null,
  serving_modes_override: null,
  colocate: null,
  serving_modes: null,
  asset_class_name: 'fasta',
  asset_group_name: 'fasta',
  genome_digest: 'g',
  names: null,
  seek_keys: null,
  ...over,
});

describe('fastaPullSize', () => {
  it("is the fasta group's default asset size", () => {
    expect(
      fastaPullSize([
        asset({ name: 'old', size: 1 }),
        asset({ name: 'new', size: 2, is_default: true }),
        asset({ asset_group_name: 'fasta_index', size: 3 }),
      ]),
    ).toBe(2);
  });

  it('falls back to the first fasta asset, and is null without one', () => {
    expect(fastaPullSize([asset({ size: 5 })])).toBe(5);
    expect(fastaPullSize([asset({ asset_group_name: 'bwa_index', size: 5 })])).toBeNull();
    expect(fastaPullSize(undefined)).toBeNull();
  });
});
