import { describe, expect, it } from 'vitest';
import { bedSnippets, downloadSnippets, regionSnippets } from './codeSnippets';
import { parseRegion } from './parseRegion';

const STORE = 'https://store.test/s/';
const DIGEST = 'EiFob05aCWgVU_B_Ae0cypnQut3cxUP1';
const SEQ = 'k3grVkjY-hoWcCUojHw6VU6GE3MZ8Sct';
const region = parseRegion('chrM:1-60', [{ name: 'chrM', length: 16569 }]);

describe('regionSnippets', () => {
  const snippets = regionSnippets({
    digest: DIGEST,
    storeUrl: STORE,
    region,
    seq: { sha512t24u: SEQ, alphabet: 'dna3bit' },
    seqFileUrl: `${STORE}sequences/k3/${SEQ}.seq`,
  });
  const code = (id: string) => snippets.find((s) => s.id === id)?.code ?? '';

  it('covers the tools with a verified remote read', () => {
    expect(snippets.map((s) => s.label)).toEqual([
      'CLI',
      'Python',
      'R',
      'Rust',
      'Node',
      'Browser JS',
    ]);
  });

  it('turns the 1-based region into 0-based, half-open coordinates', () => {
    expect(code('cli')).toBe(
      `refget store get ${SEQ} -s --start 0 --end 60 --remote ${STORE}`,
    );
    expect(code('python')).toContain(`store.get_substring("${SEQ}", 0, 60)`);
    expect(code('rust')).toContain(`store.get_substring("${SEQ}", 0, 60)?`);
    expect(code('r')).toContain('get_substring(store, seq_digest, 0, 60)');
    expect(code('r')).toContain(`get_level2(store, "${DIGEST}")`);
    expect(code('node')).toContain('encodedByteRange(0, 60, "dna3bit")');
    expect(code('browser')).toContain('store.getSubstring("chrM", 0, 60)');
  });
});

describe('downloadSnippets', () => {
  it('names the store and file, and the chosen sequences', () => {
    const all = downloadSnippets({ digest: DIGEST, storeUrl: STORE, alias: 'hg38' });
    expect(all.map((s) => s.id)).toEqual(['refget', 'python']);
    expect(all[0].code).toBe(
      `refget store export ${DIGEST} --remote ${STORE} -o hg38.fa.gz`,
    );
    const some = downloadSnippets({
      digest: DIGEST,
      storeUrl: STORE,
      alias: undefined,
      names: ['chr1', 'chrM'],
    });
    expect(some[0].code).toContain('-n chr1 -n chrM');
    expect(some[1].code).toContain(`"${DIGEST}.fa.gz", ["chr1", "chrM"])`);
  });

  it('adds refgenie only for a genome with an alias and a fasta asset', () => {
    const args = { digest: DIGEST, storeUrl: STORE, hasFastaAsset: true };
    expect(downloadSnippets({ ...args, alias: 'hg38' })[2].code).toBe(
      'refgenie init\nrefgenie pull hg38/fasta\nrefgenie seek hg38/fasta',
    );
    expect(downloadSnippets({ ...args, alias: undefined })).toHaveLength(2);
  });
});

describe('bedSnippets', () => {
  const snippets = bedSnippets({
    digest: DIGEST,
    storeUrl: STORE,
    seqTemplate: 'sequences/%s2/%s.seq',
    collectionTemplate: 'collections/%s.rgsi',
  });
  const code = (id: string) => snippets.find((s) => s.id === id)?.code ?? '';

  it('covers the same tools as one region', () => {
    expect(snippets.map((s) => s.label)).toEqual([
      'CLI',
      'Python',
      'R',
      'Rust',
      'Node',
      'Browser JS',
    ]);
  });

  it('uses each library\'s own BED export where it has one', () => {
    expect(code('cli')).toBe(
      `refget store export ${DIGEST} --bed regions.bed --remote ${STORE} -o regions.fa`,
    );
    expect(code('python')).toContain(`store.load_collection("${DIGEST}")`);
    expect(code('python')).toContain(
      `store.export_fasta_from_regions("${DIGEST}", "regions.bed", "regions.fa")`,
    );
    expect(code('r')).toContain(`get_seqs_bed_file(store, "${DIGEST}", "regions.bed", "regions.fa")`);
    expect(code('rust')).toContain(
      `store.export_fasta_from_regions("${DIGEST}", "regions.bed", "regions.fa")?;`,
    );
  });

  it('builds store URLs from the store\'s templates in the JavaScript loop', () => {
    expect(code('node')).toContain(`fetch("${STORE}collections/${DIGEST}.rgsi")`);
    expect(code('node')).toContain('fetch(`' + STORE + 'sequences/${d.slice(0, 2)}/${d}.seq`');
    expect(code('browser')).toContain(`store.openCollection("${DIGEST}")`);
    expect(code('browser')).toContain('store.getSubstring(chrom, Number(start), Number(end))');
  });
});
