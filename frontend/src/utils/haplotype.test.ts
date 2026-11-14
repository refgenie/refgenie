import { describe, expect, it } from 'vitest';
import { haplotypeOf } from './haplotype';

function of(alias: string) {
  return haplotypeOf({ aliases: [alias], synonyms: null, description: null });
}

describe('haplotypeOf', () => {
  it.each([
    // Real names from refgenie-registry stores/vgp/sources.csv.
    ['aBomVar4.paternal.1', 'paternal'],
    ['fAciRut3.2_maternal_haplotype', 'maternal'],
    ['fLatCha1.pri', 'primary'],
    ['bAcrTri1.alt', 'alternate'],
    ['bAegAlb1.hap1', 'hap1'],
    ['fArrGeo1.hap2', 'hap2'],
    ['mCamDro1.pat', 'paternal'],
    ['mCalJa1.2.pat.X', 'paternal'],
  ])('%s -> %s', (alias, expected) => {
    expect(of(alias)).toBe(expected);
  });

  it.each(['fAbrBra2.1', 'ASM4175387v1', 'mPatXyz1.1'])('%s -> null', (alias) => {
    expect(of(alias)).toBeNull();
  });

  it('falls back to synonyms, then the description', () => {
    expect(
      haplotypeOf({ aliases: ['fAbrBra2.1'], synonyms: ['x.hap2'], description: null }),
    ).toBe('hap2');
    expect(
      haplotypeOf({
        aliases: [],
        synonyms: null,
        description: 'bAcrTri1 maternal haplotype',
      }),
    ).toBe('maternal');
  });
});
