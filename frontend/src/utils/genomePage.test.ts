import { describe, expect, it } from 'vitest';
import { defaultGenomeTab } from './genomePage';

describe('defaultGenomeTab', () => {
  it('waits for the sequence list', () => {
    expect(defaultGenomeTab({ assetCount: 3, hasSequences: null })).toBeNull();
  });
  it('prefers Sequences when there are any', () => {
    expect(defaultGenomeTab({ assetCount: null, hasSequences: true })).toBe(
      'sequences',
    );
    expect(defaultGenomeTab({ assetCount: 3, hasSequences: true })).toBe('sequences');
  });
  it('falls back to Assets, then Details', () => {
    expect(defaultGenomeTab({ assetCount: null, hasSequences: false })).toBeNull();
    expect(defaultGenomeTab({ assetCount: 3, hasSequences: false })).toBe('assets');
    expect(defaultGenomeTab({ assetCount: 0, hasSequences: false })).toBe('details');
  });
});
