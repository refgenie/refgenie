/**
 * Code that does what the genome page's Download and Extract tabs do, in
 * other tools. Every call here was checked against the released packages
 * (PyPI gtars 0.11.2 and refget 0.13.0, refgenie 1.0.0a1, crates.io
 * gtars-refget 0.11.1, npm @databio/gtars 0.9.2, R gtars 0.10.1) and run
 * against hg38's chrM:1-60. Tools with no working remote-read call (the gtars
 * CLI, the @databio/gtars-node package as published) are left out.
 *
 * All of these APIs take 0-based, half-open coordinates: [start, end).
 */
import type { CodeSnippet } from '../components/common/CodeSnippets';
import type { ParsedRegion } from './parseRegion';
import { refgetExportCommand } from './fastaDownload';
import { expandTemplate } from '../services/fastaExport';

/**
 * gtars for R builds from a clone: its Rust code uses the other gtars crates
 * by path, so `remotes::install_github(subdir = "gtars-r")` cannot build it.
 */
const R_INSTALL =
  'git clone --depth 1 --branch gtars-r-v0.10.1 https://github.com/databio/gtars && R CMD INSTALL gtars/gtars-r';
const R_NOTE =
  'A shell command. Building needs Rust, and the R packages BiocGenerics, IRanges (Bioconductor) and data.table first.';

/** A string literal for Python, R, Rust or JavaScript. */
function str(value: string): string {
  return JSON.stringify(value);
}

export interface RegionSnippetArgs {
  /** The genome's sequence collection digest. */
  digest: string;
  storeUrl: string;
  region: ParsedRegion;
  /** The region's sequence: its digest and alphabet. */
  seq: { sha512t24u: string; alphabet: string };
  /** The sequence's `.seq` file URL. */
  seqFileUrl: string;
}

/** The same region, read from the genome's refget store by each tool. */
export function regionSnippets({
  storeUrl,
  digest,
  region,
  seq,
  seqFileUrl,
}: RegionSnippetArgs): CodeSnippet[] {
  const { name, start, end } = region;
  const d = seq.sha512t24u;
  const coords = `0-based, half-open: ${region.label} is [${start}, ${end}).`;
  return [
    {
      id: 'cli',
      label: 'CLI',
      install: 'pip install refget',
      code: [
        'refget store get',
        d,
        '-s',
        '--start',
        String(start),
        '--end',
        String(end),
        '--remote',
        storeUrl,
      ].join(' '),
      note: `The refget CLI, by the digest of ${name}. Coordinates are ${coords}`,
    },
    {
      id: 'python',
      label: 'Python',
      install: 'pip install gtars',
      code: [
        'from gtars.refget import RefgetStore',
        '',
        `store = RefgetStore.open_remote("refget_cache", ${str(storeUrl)})`,
        `# ${name}`,
        `print(store.get_substring(${str(d)}, ${start}, ${end}))`,
      ].join('\n'),
      note: `Reads only the bytes it needs. Coordinates are ${coords}`,
    },
    {
      id: 'r',
      label: 'R',
      install: R_INSTALL,
      code: [
        'library(gtars)',
        '',
        `store <- refget_store_open_remote("refget_cache", ${str(storeUrl)})`,
        `lvl2 <- get_level2(store, ${str(digest)})`,
        `seq_digest <- lvl2$sequences[match(${str(name)}, lvl2$names)]`,
        `get_substring(store, seq_digest, ${start}, ${end})`,
      ].join('\n'),
      note: `${R_NOTE} Coordinates are ${coords}`,
    },
    {
      id: 'rust',
      label: 'Rust',
      install: 'cargo add gtars-refget@0.11 anyhow',
      code: [
        'use gtars_refget::store::RefgetStore;',
        '',
        'fn main() -> anyhow::Result<()> {',
        `    let mut store = RefgetStore::open_remote("refget_cache", ${str(storeUrl)})?;`,
        '    store.load_sequence_index()?;',
        `    // ${name}`,
        `    println!("{}", store.get_substring(${str(d)}, ${start}, ${end})?);`,
        '    Ok(())',
        '}',
      ].join('\n'),
      note: `Coordinates are ${coords}`,
    },
    {
      id: 'node',
      label: 'Node',
      install: 'npm install @databio/gtars@0.9.2',
      code: [
        'import { readFile } from "node:fs/promises";',
        'import init, { encodedByteRange, decodeEncodedRange } from "@databio/gtars";',
        '',
        'const wasm = new URL("gtars_js_bg.wasm", import.meta.resolve("@databio/gtars"));',
        'await init({ module_or_path: await readFile(wasm) });',
        `// ${name}`,
        `const url = ${str(seqFileUrl)};`,
        `const [from, to] = encodedByteRange(${start}, ${end}, ${str(seq.alphabet)});`,
        'const res = await fetch(url, { headers: { Range: `bytes=${from}-${to - 1}` } });',
        'const bytes = new Uint8Array(await res.arrayBuffer());',
        `console.log(decodeEncodedRange(bytes, from, ${start}, ${end}, ${str(seq.alphabet)}));`,
      ].join('\n'),
      note: `The WebAssembly build, in Node 22 as an ES module (.mjs). Fetches only the bytes it needs. Coordinates are ${coords}`,
    },
    {
      id: 'browser',
      label: 'Browser JS',
      install: 'npm install @databio/gtars@0.9.2',
      code: [
        'import init, * as wasm from "@databio/gtars";',
        'import { RemoteRefgetStore } from "@databio/gtars/remote";',
        '',
        'await init();',
        `const store = new RemoteRefgetStore({ baseUrl: ${str(storeUrl)}, wasm });`,
        'await store.prepare();',
        `await store.openCollection(${str(digest)});`,
        `console.log(await store.getSubstring(${str(name)}, ${start}, ${end}));`,
      ].join('\n'),
      note: `WebAssembly in the browser. Its cache needs a Web Worker; on the main thread it still reads, without caching. Coordinates are ${coords}`,
    },
  ];
}

export interface DownloadSnippetArgs {
  digest: string;
  storeUrl: string;
  /** The genome's alias, for file names and refgenie; the digest if none. */
  alias: string | undefined;
  /** Sequences to export; all when omitted. */
  names?: string[];
  /** True if the genome has a built `fasta` asset that refgenie can pull. */
  hasFastaAsset?: boolean;
}

const DESCRIPTIONS =
  "The browser's file has names only in its headers; this one also keeps each sequence's description.";

/** The whole genome (or the chosen sequences) as FASTA, with each tool. */
export function downloadSnippets({
  digest,
  storeUrl,
  alias,
  names,
  hasFastaAsset = false,
}: DownloadSnippetArgs): CodeSnippet[] {
  const output = `${alias ?? digest}.fa.gz`;
  const snippets: CodeSnippet[] = [
    {
      id: 'refget',
      label: 'refget CLI',
      install: 'pip install refget',
      code: refgetExportCommand({ digest, storeUrl, output, names }),
      note: names
        ? 'Only the sequences chosen above. A .gz name turns on compression. ' + DESCRIPTIONS
        : 'Every sequence. A .gz name turns on compression; add -n NAME (repeatable) for only some. ' +
          DESCRIPTIONS,
    },
    {
      id: 'python',
      label: 'Python',
      install: 'pip install gtars',
      code: [
        'from gtars.refget import RefgetStore',
        '',
        `store = RefgetStore.open_remote("refget_cache", ${str(storeUrl)})`,
        names
          ? `store.export_fasta(${str(digest)}, ${str(output)}, ${pyList(names)})`
          : `store.export_fasta(${str(digest)}, ${str(output)})`,
      ].join('\n'),
      note: names
        ? 'Only the sequences chosen above.'
        : 'Every sequence; pass a list of names as a third argument for only some.',
    },
  ];
  // refgenie names the genome by alias; no alias, no snippet.
  if (hasFastaAsset && alias) {
    const target = `${alias}/fasta`;
    snippets.push({
      id: 'refgenie',
      label: 'refgenie CLI',
      install: 'pip install --pre refgenie',
      code: ['refgenie init', `refgenie pull ${target}`, `refgenie seek ${target}`].join(
        '\n',
      ),
      note: 'Pulls the built fasta asset (.fa with its .fai index) from a refgenie server; seek prints where it is.',
    });
  }
  return snippets;
}

function pyList(names: string[]): string {
  return `[${names.map(str).join(', ')}]`;
}

export interface BedSnippetArgs {
  /** The genome's sequence collection digest. */
  digest: string;
  storeUrl: string;
  /** The store's `.seq` path template, e.g. `sequences/%s2/%s.seq`. */
  seqTemplate: string;
  /** The store's collection index path template, e.g. `collections/%s.rgsi`. */
  collectionTemplate: string;
}

/** A store path template as a JavaScript template-literal body, over digest `d`. */
function templateLiteral(template: string): string {
  return template
    .split('%s2')
    .join('${d.slice(0, 2)}')
    .split('%s4')
    .join('${d.slice(0, 4)}')
    .split('%s')
    .join('${d}');
}

const NATIVE_BED_NOTE =
  'Writes regions.fa with chr:start-end headers (BED coordinates); it skips the name column. ' +
  'The BED must be tab-separated, with no track, browser or # lines; a .bed.gz works too.';

/**
 * A BED file of regions, read from the genome's refget store into a FASTA
 * file, with each tool. Every call here was run against hg38 with the
 * released packages named in this file's header.
 */
export function bedSnippets({
  digest,
  storeUrl,
  seqTemplate,
  collectionTemplate,
}: BedSnippetArgs): CodeSnippet[] {
  const base = storeUrl.endsWith('/') ? storeUrl : `${storeUrl}/`;
  const collectionUrl = new URL(expandTemplate(collectionTemplate, digest), base).href;
  return [
    {
      id: 'cli',
      label: 'CLI',
      install: 'pip install refget',
      code: `refget store export ${digest} --bed regions.bed --remote ${storeUrl} -o regions.fa`,
      note: NATIVE_BED_NOTE,
    },
    {
      id: 'python',
      label: 'Python',
      install: 'pip install gtars',
      code: [
        'from gtars.refget import RefgetStore',
        '',
        `store = RefgetStore.open_remote("refget_cache", ${str(storeUrl)})`,
        `store.load_collection(${str(digest)})`,
        `store.export_fasta_from_regions(${str(digest)}, "regions.bed", "regions.fa")`,
      ].join('\n'),
      note: `Reads only the bytes it needs. ${NATIVE_BED_NOTE}`,
    },
    {
      id: 'r',
      label: 'R',
      install: R_INSTALL,
      code: [
        'library(gtars)',
        '',
        `store <- refget_store_open_remote("refget_cache", ${str(storeUrl)})`,
        `invisible(get_collection(store, ${str(digest)}))`,
        `get_seqs_bed_file(store, ${str(digest)}, "regions.bed", "regions.fa")`,
      ].join('\n'),
      note: `${R_NOTE} ${NATIVE_BED_NOTE}`,
    },
    {
      id: 'rust',
      label: 'Rust',
      install: 'cargo add gtars-refget@0.11 anyhow',
      code: [
        'use gtars_refget::store::RefgetStore;',
        '',
        'fn main() -> anyhow::Result<()> {',
        `    let mut store = RefgetStore::open_remote("refget_cache", ${str(storeUrl)})?;`,
        `    store.load_collection(${str(digest)})?;`,
        `    store.export_fasta_from_regions(${str(digest)}, "regions.bed", "regions.fa")?;`,
        '    Ok(())',
        '}',
      ].join('\n'),
      note: NATIVE_BED_NOTE,
    },
    {
      id: 'node',
      label: 'Node',
      install: 'npm install @databio/gtars@0.9.2',
      code: [
        'import { readFile, writeFile } from "node:fs/promises";',
        'import init, { encodedByteRange, decodeEncodedRange } from "@databio/gtars";',
        'import { parseRgsi } from "@databio/gtars/remote";',
        '',
        'const wasm = new URL("gtars_js_bg.wasm", import.meta.resolve("@databio/gtars"));',
        'await init({ module_or_path: await readFile(wasm) });',
        '',
        `const index = await (await fetch(${str(collectionUrl)})).text();`,
        'const seqs = new Map(parseRgsi(index).map((s) => [s.name, s]));',
        '',
        'let fasta = "";',
        'for (const line of (await readFile("regions.bed", "utf8")).split("\\n")) {',
        '  if (!line.trim() || /^(#|track|browser)/.test(line)) continue;',
        '  const [chrom, from, to, name] = line.trim().split(/\\s+/);',
        '  const { sha512t24u: d, alphabet } = seqs.get(chrom);',
        '  const [start, end] = [Number(from), Number(to)];',
        '  const [bs, be] = encodedByteRange(start, end, alphabet);',
        `  const res = await fetch(\`${base}${templateLiteral(seqTemplate)}\`, {`,
        '    headers: { Range: `bytes=${bs}-${be - 1}` },',
        '  });',
        '  const bytes = new Uint8Array(await res.arrayBuffer());',
        '  const coords = `${chrom}:${start}-${end}`;',
        '  fasta += `>${name ? `${name}::${coords}` : coords}\\n`;',
        '  fasta += `${decodeEncodedRange(bytes, bs, start, end, alphabet)}\\n`;',
        '}',
        'await writeFile("regions.fa", fasta);',
      ].join('\n'),
      note: 'The WebAssembly build, in Node 22 as an ES module (.mjs). Fetches only the bytes it needs, one region at a time. Headers match this page: name::chr:start-end.',
    },
    {
      id: 'browser',
      label: 'Browser JS',
      install: 'npm install @databio/gtars@0.9.2',
      code: [
        'import init, * as wasm from "@databio/gtars";',
        'import { RemoteRefgetStore } from "@databio/gtars/remote";',
        '',
        'await init();',
        `const store = new RemoteRefgetStore({ baseUrl: ${str(storeUrl)}, wasm });`,
        'await store.prepare();',
        `await store.openCollection(${str(digest)});`,
        '',
        '// file: a BED File, e.g. from <input type="file">',
        'async function bedToFasta(file) {',
        '  let fasta = "";',
        '  for (const line of (await file.text()).split("\\n")) {',
        '    if (!line.trim() || /^(#|track|browser)/.test(line)) continue;',
        '    const [chrom, start, end, name] = line.trim().split(/\\s+/);',
        '    const bases = await store.getSubstring(chrom, Number(start), Number(end));',
        '    const coords = `${chrom}:${start}-${end}`;',
        '    fasta += `>${name ? `${name}::${coords}` : coords}\\n${bases}\\n`;',
        '  }',
        '  return fasta;',
        '}',
      ].join('\n'),
      note: 'WebAssembly in the browser; returns the FASTA text. Its cache needs a Web Worker; on the main thread it still reads, without caching. Headers match this page.',
    },
  ];
}
