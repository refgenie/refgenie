import { useEffect, useMemo, useRef, useState } from 'react';
import type { DragEvent } from 'react';
import { DigestChip } from '../common/DigestChip';
import { ScopedLink } from '../common/ScopedLink';
import { ComparisonResult } from './ComparisonResult';
import { useResource } from '../../hooks/useResource';
import { useApiClient } from '../../hooks/useApiClient';
import { qk } from '../../services/queryKeys';
import { NO_RETRY } from '../../services/resourceCache';
import { getGenome } from '../../services/resources/genomes';
import { FASTA_ACCEPT, looksLikeFasta } from '../../vendor/fastaDigest/fastaDigest.ts';
import type { FastaDigestResult } from '../../vendor/fastaDigest/fastaDigest.ts';
import { digestFile, isAbortError, type WorkerFactory } from '../../vendor/fastaDigest/fastaDigestClient.ts';
import {
  comparisonArrays,
  compareSeqcols,
  summarizeFileComparison,
} from '../../utils/seqcolLocal';
import { formatBytes, formatDigest } from '../../utils/format';
import { preferredAlias } from '../../utils/genomes';
import { listAliases } from '../../services/resources/aliases';
import { cn } from '../../utils/cn';
import type { SequenceRow } from './types';

type State =
  | { kind: 'idle' }
  | { kind: 'reading'; file: File; bytesRead: number }
  | { kind: 'done'; file: File; result: FastaDigestResult; seconds: number }
  | { kind: 'error'; file: File | null; message: string };

export interface CompareFastaProps {
  /** This genome's digest. */
  digest: string;
  /** This genome's sequences, in collection order; undefined while loading. */
  sequences: SequenceRow[] | undefined;
  /** Tests pass a fake worker. */
  workerFactory?: WorkerFactory;
}

/** Another genome in refgenie whose digest equals the file's, and its name. */
function useMatchingGenome(fileDigest: string | undefined): {
  digest: string;
  label: string;
} | null {
  const client = useApiClient();
  const genome = useResource(
    qk.genome(fileDigest ?? ''),
    ({ signal }) => getGenome(client, fileDigest as string, { signal }),
    { enabled: !!fileDigest, retry: NO_RETRY },
  );
  const found = genome.data?.digest;
  const params = { genome_digest: found, pageSize: 100 };
  const aliases = useResource(
    qk.aliases(params),
    ({ signal }) => listAliases(client, params, { signal }),
    { enabled: !!found },
  );
  if (!found || !aliases.data) return null;
  const label = preferredAlias(aliases.data.results.map((alias) => alias.name));
  return { digest: found, label: label ?? formatDigest(found) };
}

/**
 * The Compare tab's "your own FASTA" box: drop or pick a FASTA file, digest it
 * in the browser (a worker running the gtars wasm hasher), and compare its
 * sequence collection with this genome's.
 */
export function CompareFasta({ digest, sequences, workerFactory }: CompareFastaProps) {
  const [state, setState] = useState<State>({ kind: 'idle' });
  const [dragging, setDragging] = useState(false);
  const cancelRef = useRef<(() => void) | null>(null);
  const inputRef = useRef<HTMLInputElement | null>(null);

  // Stop a running digest when the tab closes. Clearing the ref first makes
  // the job's handlers skip their setState on the unmounted component.
  useEffect(
    () => () => {
      const stop = cancelRef.current;
      cancelRef.current = null;
      stop?.();
    },
    [],
  );

  const start = (file: File) => {
    cancelRef.current?.();
    if (!looksLikeFasta(file.name) && !/\.b?gz$/i.test(file.name)) {
      setState({
        kind: 'error',
        file,
        message: 'Choose a FASTA file: .fa, .fasta or .fna, optionally gzipped (.gz).',
      });
      return;
    }
    setState({ kind: 'reading', file, bytesRead: 0 });
    const began = performance.now();
    const job = digestFile(
      file,
      {
        onProgress: (bytesRead) =>
          cancelRef.current === job.cancel &&
          setState((s) =>
            s.kind === 'reading' && s.file === file ? { ...s, bytesRead } : s,
          ),
      },
      workerFactory,
    );
    cancelRef.current = job.cancel;
    job.done.then(
      ({ result }) => {
        if (cancelRef.current !== job.cancel) return;
        cancelRef.current = null;
        setState({
          kind: 'done',
          file,
          result,
          seconds: (performance.now() - began) / 1000,
        });
      },
      (err: unknown) => {
        if (cancelRef.current !== job.cancel) return;
        cancelRef.current = null;
        if (isAbortError(err)) {
          setState({ kind: 'idle' });
          return;
        }
        setState({
          kind: 'error',
          file,
          message: err instanceof Error ? err.message : String(err),
        });
      },
    );
  };

  const cancel = () => {
    const stop = cancelRef.current;
    cancelRef.current = null;
    stop?.();
    setState({ kind: 'idle' });
  };

  const onDrop = (event: DragEvent<HTMLDivElement>) => {
    event.preventDefault();
    setDragging(false);
    const file = event.dataTransfer.files?.[0];
    if (file && state.kind !== 'reading') start(file);
  };

  const done = state.kind === 'done' ? state : null;
  const fileRows = done?.result.sequences;
  const comparison = useMemo(
    () =>
      done && sequences
        ? compareSeqcols(
            comparisonArrays(sequences),
            comparisonArrays(done.result.sequences),
            {
              a: digest,
              b: done.result.digest,
            },
          )
        : null,
    [done, sequences, digest],
  );
  const summary = useMemo(
    () =>
      done && sequences && fileRows
        ? summarizeFileComparison(sequences, fileRows, {
            genome: digest,
            file: done.result.digest,
          })
        : null,
    [done, sequences, fileRows, digest],
  );
  const lookup = done && done.result.digest !== digest ? done.result.digest : undefined;
  const match = useMatchingGenome(lookup);

  const reading = state.kind === 'reading' ? state : null;
  const percent =
    reading && reading.file.size > 0
      ? Math.min(100, Math.round((100 * reading.bytesRead) / reading.file.size))
      : 0;

  return (
    <section className="rg-flat-section">
      <h2>Compare your own FASTA</h2>
      <div className="flex flex-col gap-4">
        <p className="rg-muted text-sm">
          See how a FASTA file on your computer relates to this genome. Your browser
          reads the file and computes its sequence collection digest; the file never
          leaves your computer.
        </p>
        <div
          className={cn('rg-dropzone', dragging && 'rg-dropzone--active')}
          onDragOver={(event) => {
            event.preventDefault();
            if (!reading) setDragging(true);
          }}
          onDragLeave={() => setDragging(false)}
          onDrop={onDrop}
          data-testid="fasta-dropzone"
        >
          {reading ? (
            <div className="rg-dropzone__body">
              <p className="text-sm">
                Reading{' '}
                <strong className="rg-dropzone__name">{reading.file.name}</strong>
              </p>
              <progress
                className="rg-progress__bar"
                value={percent}
                max={100}
                aria-label={`Reading ${reading.file.name}: ${percent}%`}
              />
              <p className="rg-muted text-sm">
                {formatBytes(reading.bytesRead)} of {formatBytes(reading.file.size)} (
                {percent}
                %)
              </p>
              <button type="button" className="rg-btn rg-btn--sm" onClick={cancel}>
                Cancel
              </button>
            </div>
          ) : (
            <div className="rg-dropzone__body">
              <p>Drop a FASTA file here</p>
              <p className="rg-muted text-sm">.fa, .fasta or .fna, plain or gzipped</p>
              <button
                type="button"
                className="rg-btn rg-btn--sm"
                onClick={() => inputRef.current?.click()}
              >
                Choose a file…
              </button>
            </div>
          )}
          <input
            ref={inputRef}
            type="file"
            accept={FASTA_ACCEPT}
            className="sr-only"
            aria-label="FASTA file to compare"
            tabIndex={-1}
            onChange={(event) => {
              const file = event.target.files?.[0];
              event.target.value = '';
              if (file) start(file);
            }}
          />
        </div>

        {state.kind === 'error' && (
          <p className="rg-state rg-state--error" role="alert">
            {state.file ? `${state.file.name}: ` : ''}
            {state.message}
          </p>
        )}

        {done && summary && comparison && (
          <div className="flex flex-col gap-4" aria-live="polite">
            <div className="flex flex-col gap-2">
              <h3 className="text-base font-semibold">{summary.headline}</h3>
              <ul className="rg-compare-fasta__details flex flex-col gap-1">
                {summary.details.map((line) => (
                  <li key={line}>{line}</li>
                ))}
              </ul>
              {match && (
                <p>
                  Your file matches{' '}
                  <ScopedLink
                    className="rg-link"
                    to={`/genomes/${encodeURIComponent(match.digest)}`}
                  >
                    {match.label}
                  </ScopedLink>{' '}
                  in refgenie.
                </p>
              )}
              <p className="rg-muted text-sm flex items-center gap-2 flex-wrap">
                <span>
                  {done.file.name} · {formatBytes(done.file.size)} ·{' '}
                  {done.result.n_sequences.toLocaleString()}{' '}
                  {done.result.n_sequences === 1 ? 'sequence' : 'sequences'} · read in{' '}
                  {done.seconds < 10
                    ? done.seconds.toFixed(1)
                    : Math.round(done.seconds)}{' '}
                  s
                </span>
              </p>
              <p className="text-sm flex items-center gap-2 flex-wrap">
                <span>Your file&apos;s digest:</span>
                <DigestChip digest={done.result.digest} length={32} />
              </p>
            </div>
            <ComparisonResult
              comparison={comparison}
              caption={<>Compared with your file</>}
              otherName="your file"
              otherColumn="In your file"
            />
          </div>
        )}
      </div>
    </section>
  );
}
