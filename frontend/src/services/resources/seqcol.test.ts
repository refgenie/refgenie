import { describe, expect, it, vi } from 'vitest';
import { coerceSeqcolServices } from '../config';
import { findCollectionsSharing, seqcolBase, serviceForStore } from './seqcol';

const SERVICES = [
  {
    store: 'jungle',
    url: '/seqcol',
    refget_store_url: 'https://s3/refget-store/jungle/',
  },
  {
    store: 'vgp',
    url: '/seqcol/vgp',
    refget_store_url: 'https://s3/refget-store/vgp/',
  },
];

describe('coerceSeqcolServices', () => {
  it('reads seqcol.services from service-info, dropping malformed entries', () => {
    expect(
      coerceSeqcolServices({
        services: [
          {
            store: 'vgp',
            url: '/seqcol/vgp',
            refget_store: { url: 'https://s3/vgp/' },
          },
          { url: 42 },
        ],
      }),
    ).toEqual([
      { store: 'vgp', url: '/seqcol/vgp', refget_store_url: 'https://s3/vgp/' },
    ]);
    expect(coerceSeqcolServices(undefined)).toEqual([]);
  });
});

describe('serviceForStore', () => {
  it("picks the service that reads the genome's store, else the first", () => {
    expect(serviceForStore(SERVICES, 'https://s3/refget-store/vgp')?.store).toBe('vgp');
    expect(serviceForStore(SERVICES, null)?.store).toBe('jungle');
  });

  it('builds the base on the API server', () => {
    expect(seqcolBase('https://api.test/v4', SERVICES[1])).toBe(
      'https://api.test/seqcol/vgp',
    );
  });
});

describe('findCollectionsSharing', () => {
  it('asks every store and merges the digests, skipping stores that fail', async () => {
    const fetchMock = vi.fn(async (url: string) => {
      if (url.includes('/seqcol/vgp/')) return new Response('nope', { status: 500 });
      return new Response(JSON.stringify({ results: ['a', 'b'] }), {
        headers: { 'content-type': 'application/json' },
      });
    });
    vi.stubGlobal('fetch', fetchMock);
    try {
      const found = await findCollectionsSharing(
        'https://api.test/v4',
        SERVICES,
        'lengths',
        'L',
      );
      expect(found).toEqual(['a', 'b']);
      expect(String(fetchMock.mock.calls[0][0])).toBe(
        'https://api.test/seqcol/list/collection?lengths=L&page_size=100',
      );
    } finally {
      vi.unstubAllGlobals();
    }
  });
});
