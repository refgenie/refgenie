import { describe, expect, it } from 'vitest';
import { accessionUrl, safeHttpUrl, taxonUrl } from './links';

describe('safeHttpUrl', () => {
  it.each(['https://example.org/a', 'http://example.org', '/openapi.json', 'logo.png'])(
    'keeps %j',
    (url) => expect(safeHttpUrl(url)).toBe(url),
  );

  it.each([
    'javascript:alert(1)',
    '  javascript:alert(1)',
    'JAVASCRIPT:alert(1)',
    'java\tscript:alert(1)',
    'data:text/html,<script>alert(1)</script>',
    'vbscript:msgbox',
    'file:///etc/passwd',
    '',
    null,
    undefined,
  ])('rejects %j', (url) => expect(safeHttpUrl(url)).toBeUndefined());
});

describe('accessionUrl', () => {
  it('uses a safe FHR url, else falls back to NCBI', () => {
    expect(accessionUrl('GCA_1', 'https://x.org/a')).toBe('https://x.org/a');
    expect(accessionUrl('GCA_1', 'javascript:alert(1)')).toBe(
      'https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_1/',
    );
    expect(accessionUrl('other', 'javascript:alert(1)')).toBeUndefined();
  });
});

describe('taxonUrl', () => {
  it('drops an unsafe taxon_uri', () => {
    expect(taxonUrl('javascript:alert(1)', 9606)).toBe('https://identifiers.org/taxonomy:9606');
    expect(taxonUrl('javascript:alert(1)', null)).toBeUndefined();
  });
});
