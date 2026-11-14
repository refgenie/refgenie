import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import {
  assertDecodable,
  fastaChunks,
  fastaSeqsFromRows,
  fastaSize,
  isDecodable,
  loadCollection,
  loadGtars,
  type FastaCollection,
} from '../services/fastaExport';
import { saveStream } from '../services/saveStream';
import { bgzf } from '../services/bgzf';
import { useGenomeSequences } from './queries/useGenomeSequences';
import {
  estimateFastaGzBytes,
  primaryChromosomes,
} from '../utils/fastaDownload';

export type FastaDownloadStatus =
  | { kind: 'idle' }
  | { kind: 'saving'; fraction: number; written: number }
  | { kind: 'saved' }
  | { kind: 'error'; message: string };

export interface FastaDownloadState {
  /** The full sequence list, once loaded. */
  collection: FastaCollection | null;
  status: FastaDownloadStatus;
  /** False: all sequences. True: only `selected`. */
  choosing: boolean;
  setChoosing: (choosing: boolean) => void;
  selected: Set<string>;
  setSelected: (selected: Set<string>) => void;
  /** Sequences the browser cannot rebuild, each with a short reason. */
  unavailable: ReadonlyMap<string, string>;
  /** Sequences that will be written, or null while loading. */
  count: number | null;
  /** Uncompressed FASTA size of the selection. */
  totalBytes: number | null;
  /** Likely .fa.gz size of the selection (an estimate). */
  estimatedBytes: number | null;
  /** The chosen sequence names, for the matching commands; undefined for all. */
  names: string[] | undefined;
  start: () => Promise<void>;
  cancel: () => void;
}

/**
 * The genome page's FASTA download: loads the sequence list, holds the
 * sequence choice, and runs the download. It lives on the page, not in the
 * dialog, so closing the dialog does not stop a download.
 */
export function useFastaDownload({
  storeUrl,
  digest,
  fileStem,
}: {
  /** Null when the genome has no refget store: nothing loads. */
  storeUrl: string | null;
  digest: string;
  fileStem: string;
}): FastaDownloadState {
  const [collection, setCollection] = useState<FastaCollection | null>(null);
  const [status, setStatus] = useState<FastaDownloadStatus>({ kind: 'idle' });
  const [choosing, setChoosing] = useState(false);
  const [chosen, setChosen] = useState<Set<string> | null>(null);
  const controller = useRef<AbortController | null>(null);

  // The Sequences tab loads the same list; reuse it and skip the .rgsi.
  const sequences = useGenomeSequences(storeUrl ? digest : undefined);
  const known = useMemo(
    () => (sequences.data ? fastaSeqsFromRows(sequences.data.sequences) : null),
    [sequences.data],
  );
  const waiting = sequences.isPending;

  // The manifest, plus the .rgsi if the API had no list: small, and gives
  // the exact file size up front.
  useEffect(() => {
    if (!storeUrl || waiting) return;
    let stale = false;
    loadCollection(storeUrl, digest, undefined, known)
      .then((loaded) => {
        if (!stale) setCollection(loaded);
      })
      .catch((error: unknown) => {
        if (!stale) setStatus({ kind: 'error', message: messageOf(error) });
      });
    return () => {
      stale = true;
    };
  }, [storeUrl, digest, known, waiting]);

  // Leaving the page stops a download in progress.
  useEffect(() => () => controller.current?.abort(), []);

  // Sequences the browser cannot rebuild: disabled in the picker, and only
  // a chosen one stops start(), so they never block the rest.
  const unavailable = useMemo(() => {
    const map = new Map<string, string>();
    for (const seq of collection?.seqs ?? []) {
      if (!isDecodable(seq)) {
        map.set(seq.name, `${seq.alphabet || 'unknown'} alphabet, not supported here`);
      }
    }
    return map;
  }, [collection]);

  // Starts with the primary chromosomes when the names have them, else everything.
  const selected = useMemo(() => {
    if (chosen || !collection) return chosen ?? new Set<string>();
    const names = collection.seqs
      .filter((seq) => !unavailable.has(seq.name))
      .map((seq) => seq.name);
    return new Set(primaryChromosomes(names) ?? names);
  }, [chosen, collection, unavailable]);

  const exported = useMemo(() => {
    if (!collection) return null;
    if (!choosing) return collection;
    return {
      ...collection,
      seqs: collection.seqs.filter(
        (seq) => selected.has(seq.name) && !unavailable.has(seq.name),
      ),
    };
  }, [collection, choosing, selected, unavailable]);

  const totalBytes = exported ? fastaSize(exported.seqs) : null;
  const estimatedBytes = exported
    ? estimateFastaGzBytes(exported.seqs.reduce((sum, seq) => sum + seq.length, 0))
    : null;
  const output = `${fileStem}.fa.gz`;
  const names = useMemo(
    () => (choosing && exported ? exported.seqs.map((seq) => seq.name) : undefined),
    [choosing, exported],
  );

  const start = useCallback(async () => {
    if (!storeUrl || !exported || totalBytes === null || exported.seqs.length === 0)
      return;
    if (controller.current) return;
    try {
      assertDecodable(exported.seqs);
    } catch (error) {
      setStatus({ kind: 'error', message: messageOf(error) });
      return;
    }
    const abort = new AbortController();
    controller.current = abort;
    let fraction = 0;
    let written = 0;
    setStatus({ kind: 'saving', fraction, written });
    try {
      // saveStream's first await opens the save picker, inside the click gesture.
      await saveStream({
        suggestedName: output,
        // Uncompressed, so an over-estimate: only rules out a too-large Blob.
        estimatedBytes: totalBytes,
        fileType: {
          description: 'Compressed FASTA',
          accept: { 'application/gzip': ['.gz'] },
        },
        signal: abort.signal,
        chunks: async function* () {
          const wasm = await loadGtars();
          const fasta = fastaChunks({
            baseUrl: storeUrl,
            digest,
            wasm,
            collection: exported,
            signal: abort.signal,
            onProgress: (done, total) => {
              fraction = total ? done / total : 1;
            },
          });
          let shown = 0;
          for await (const chunk of bgzf(fasta)) {
            written += chunk.length;
            // At most a few renders a second, however small the chunks.
            if (performance.now() - shown > 250) {
              shown = performance.now();
              setStatus({ kind: 'saving', fraction, written });
            }
            yield chunk;
          }
        },
      });
      setStatus({ kind: 'saved' });
    } catch (error) {
      if (isAbort(error)) setStatus({ kind: 'idle' });
      else setStatus({ kind: 'error', message: messageOf(error) });
    } finally {
      controller.current = null;
    }
  }, [storeUrl, digest, exported, totalBytes, output]);

  const cancel = useCallback(() => controller.current?.abort(), []);

  return {
    collection,
    status,
    choosing,
    setChoosing,
    selected,
    setSelected: setChosen,
    unavailable,
    count: exported ? exported.seqs.length : null,
    totalBytes,
    estimatedBytes,
    names,
    start,
    cancel,
  };
}

function isAbort(error: unknown): boolean {
  return error instanceof DOMException && error.name === 'AbortError';
}

function messageOf(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}
