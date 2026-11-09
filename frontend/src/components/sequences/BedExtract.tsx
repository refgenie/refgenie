import { useEffect, useId, useMemo, useRef, useState } from 'react';
import type { DragEvent } from 'react';
import { Icon } from '../common/Icon';
import { Copyable } from '../common/Copyable';
import { cn } from '../../utils/cn';
import { BED_CONCURRENCY, extractBed } from '../../services/bedExtract';
import { isExtractable, type RegionSeq } from '../../services/regionExtract';
import {
  MAX_BED_BASES,
  MAX_BED_REGIONS,
  bedFastaFileName,
  bedOverCap,
  exampleBed,
  parseBed,
  readBedFile,
} from '../../utils/parseBed';
import type { SequenceRow } from './types';

/** How much of the FASTA to show on the page; copy and download get it all. */
const PREVIEW_CHARS = 20_000;
/** How many problem lines to list. */
const SHOWN_PROBLEMS = 20;

export interface BedExtractProps {
  storeUrl: string;
  sequences: SequenceRow[];
}

type Source = { kind: 'paste' } | { kind: 'file'; name: string; text: string };

type Run =
  | { kind: 'idle' }
  | { kind: 'working'; done: number; total: number }
  | { kind: 'cancelled' }
  | { kind: 'error'; message: string }
  | { kind: 'done'; fasta: string; records: number; bases: number; fileName: string };

const fmt = (n: number) => n.toLocaleString();
const readable = (seq: { alphabet?: string | null }) => isExtractable(seq.alphabet);

/**
 * The Extract tab's BED mode: drop, choose or paste a BED file, see what it
 * holds, then read every region in the browser into one FASTA file.
 */
export function BedExtract({ storeUrl, sequences }: BedExtractProps) {
  const pasteId = useId();
  const hintId = useId();
  const fileInput = useRef<HTMLInputElement | null>(null);
  const [source, setSource] = useState<Source>({ kind: 'paste' });
  const [pasted, setPasted] = useState('');
  const [dragging, setDragging] = useState(false);
  const [readError, setReadError] = useState<string | null>(null);
  const [run, setRun] = useState<Run>({ kind: 'idle' });
  const controller = useRef<AbortController | null>(null);

  const text = source.kind === 'file' ? source.text : pasted;
  const bed = useMemo(() => parseBed(text, sequences, readable), [text, sequences]);
  const overCap = bedOverCap(bed);
  const example = useMemo(() => exampleBed(sequences.filter(readable)), [sequences]);
  const working = run.kind === 'working';

  useEffect(() => () => controller.current?.abort(), []);

  const reset = () => {
    controller.current?.abort();
    setRun({ kind: 'idle' });
    setReadError(null);
  };

  const loadFile = async (file: File) => {
    reset();
    try {
      const fileText = await readBedFile(file);
      setSource({ kind: 'file', name: file.name, text: fileText });
    } catch (error) {
      setReadError(error instanceof Error ? error.message : String(error));
    }
  };

  const onDrop = (event: DragEvent<HTMLDivElement>) => {
    event.preventDefault();
    setDragging(false);
    const file = event.dataTransfer.files?.[0];
    if (file) void loadFile(file);
  };

  const start = async () => {
    controller.current?.abort();
    const abort = new AbortController();
    controller.current = abort;
    const byName = new Map<string, RegionSeq>();
    for (const row of sequences) {
      byName.set(row.name, {
        name: row.name,
        length: row.length,
        alphabet: row.alphabet ?? null,
        sha512t24u: row.sha512t24u,
      });
    }
    const total = bed.regions.length;
    setRun({ kind: 'working', done: 0, total });
    try {
      const records = await extractBed({
        baseUrl: storeUrl,
        regions: bed.regions,
        sequences: byName,
        signal: abort.signal,
        onProgress: (done) => {
          if (!abort.signal.aborted) setRun({ kind: 'working', done, total });
        },
      });
      if (abort.signal.aborted) return;
      setRun({
        kind: 'done',
        fasta: records.join(''),
        records: records.length,
        bases: bed.bases,
        fileName: bedFastaFileName(source.kind === 'file' ? source.name : null),
      });
    } catch (error) {
      if (abort.signal.aborted) return;
      setRun({
        kind: 'error',
        message: error instanceof Error ? error.message : String(error),
      });
    }
  };

  const cancel = () => {
    controller.current?.abort();
    setRun({ kind: 'cancelled' });
  };

  const download = () => {
    if (run.kind !== 'done') return;
    const url = URL.createObjectURL(new Blob([run.fasta], { type: 'text/plain' }));
    const link = document.createElement('a');
    link.href = url;
    link.download = run.fileName;
    document.body.appendChild(link);
    link.click();
    link.remove();
    setTimeout(() => URL.revokeObjectURL(url), 0);
  };

  const preview = run.kind === 'done' ? previewOf(run.fasta) : null;
  const hasInput = text.trim() !== '';

  return (
    <section className="rg-bed" aria-label="Extract regions from a BED file">
      <div
        className={cn('rg-bed__drop', dragging && 'rg-bed__drop--over')}
        data-testid="bed-drop"
        onDragOver={(event) => {
          event.preventDefault();
          if (!working) setDragging(true);
        }}
        onDragLeave={() => setDragging(false)}
        onDrop={(event) => {
          if (working) {
            event.preventDefault();
            return;
          }
          onDrop(event);
        }}
      >
        {source.kind === 'file' ? (
          <div className="rg-bed__file">
            <Icon name="file" />
            <span className="rg-bed__file-name">{source.name}</span>
            <button
              type="button"
              className="rg-btn rg-btn--sm"
              disabled={working}
              onClick={() => {
                reset();
                setSource({ kind: 'paste' });
                if (fileInput.current) fileInput.current.value = '';
              }}
            >
              Remove
            </button>
          </div>
        ) : (
          <>
            <p className="rg-bed__drop-text">
              Drop a BED file here (.bed or .bed.gz), or{' '}
              <button
                type="button"
                className="rg-btn rg-btn--bare"
                disabled={working}
                onClick={() => fileInput.current?.click()}
              >
                choose a file
              </button>
              . Or paste BED lines:
            </p>
            <label className="sr-only" htmlFor={pasteId}>
              BED lines
            </label>
            <textarea
              id={pasteId}
              className="rg-field__input rg-field__input--mono rg-bed__paste"
              rows={4}
              spellCheck={false}
              placeholder={example ?? 'chr1\t999999\t1000500\tname'}
              aria-describedby={hintId}
              value={pasted}
              disabled={working}
              onChange={(event) => {
                reset();
                setPasted(event.target.value);
              }}
            />
          </>
        )}
        <input
          ref={fileInput}
          type="file"
          className="sr-only"
          tabIndex={-1}
          aria-label="Choose a BED file"
          accept=".bed,.gz,.txt,.tsv,text/plain"
          onChange={(event) => {
            const file = event.target.files?.[0];
            if (file) void loadFile(file);
          }}
        />
      </div>
      <p id={hintId} className="rg-region__hint">
        BED coordinates: 0-based, end not included. Each region becomes one FASTA record
        named like bedtools getfasta -name: <code>name::chr:start-end</code>, or{' '}
        <code>chr:start-end</code> without a name column. Strand is ignored. Up to{' '}
        {fmt(MAX_BED_REGIONS)} regions and {fmt(MAX_BED_BASES)} bases in the browser;
        for more, use the code below.
      </p>

      {readError && (
        <p className="rg-region__error" role="alert">
          {readError}
        </p>
      )}

      {hasInput && (
        <div className="rg-bed__summary" aria-live="polite">
          <p className="rg-bed__count">
            <strong>
              {fmt(bed.regions.length)}{' '}
              {bed.regions.length === 1 ? 'region' : 'regions'}
            </strong>
            {' · '}
            {fmt(bed.bases)} {bed.bases === 1 ? 'base' : 'bases'}
            {bed.problems.length > 0 && (
              <>
                {' · '}
                <span className="rg-bed__bad">
                  {fmt(bed.problems.length)}{' '}
                  {bed.problems.length === 1 ? 'line' : 'lines'} with problems, left out
                </span>
              </>
            )}
            {bed.skipped > 0 && (
              <span className="rg-muted">
                {' · '}
                {fmt(bed.skipped)} header or blank{' '}
                {bed.skipped === 1 ? 'line' : 'lines'} skipped
              </span>
            )}
          </p>
          {bed.problems.length > 0 && (
            <ul className="rg-bed__problems" aria-label="Problem lines">
              {bed.problems.slice(0, SHOWN_PROBLEMS).map((p) => (
                <li key={p.line}>
                  <span className="rg-bed__line">Line {p.line}</span>{' '}
                  <code className="rg-bed__line-text">{shorten(p.text)}</code>{' '}
                  {p.message}
                </li>
              ))}
              {bed.problems.length > SHOWN_PROBLEMS && (
                <li className="rg-muted">
                  and {fmt(bed.problems.length - SHOWN_PROBLEMS)} more
                </li>
              )}
            </ul>
          )}
          {overCap && (
            <p className="rg-region__error" role="alert">
              {overCap} That is too much to read in the browser; use one of the tools
              below, or split the file.
            </p>
          )}
        </div>
      )}

      <div className="rg-bed__actions">
        {working ? (
          <>
            <progress
              className="rg-bed__progress"
              value={run.done}
              max={run.total}
              aria-label="Regions read"
            />
            <span className="rg-region__meta">
              {fmt(run.done)} of {fmt(run.total)} regions
            </span>
            <button type="button" className="rg-btn" onClick={cancel}>
              Cancel
            </button>
          </>
        ) : (
          <>
            <button
              type="button"
              className="rg-btn rg-btn--primary"
              disabled={bed.regions.length === 0 || !!overCap}
              onClick={() => void start()}
            >
              {bed.regions.length > 0
                ? `Extract ${fmt(bed.regions.length)} ${bed.regions.length === 1 ? 'region' : 'regions'}`
                : 'Extract'}
            </button>
            {example && (
              <button
                type="button"
                className="rg-btn"
                onClick={() => {
                  reset();
                  setSource({ kind: 'paste' });
                  setPasted(example);
                }}
              >
                Try an example
              </button>
            )}
          </>
        )}
      </div>
      {working && (
        <p className="rg-region__hint">
          Reading {BED_CONCURRENCY} regions at a time; only their bytes are downloaded.
        </p>
      )}

      {run.kind === 'cancelled' && <p className="rg-region__hint">Cancelled.</p>}
      {run.kind === 'error' && (
        <p className="rg-region__error" role="alert">
          {run.message}
        </p>
      )}
      {run.kind === 'done' && preview && (
        <div className="rg-region__out">
          <div className="rg-region__bar">
            <span className="rg-region__meta">
              {fmt(run.records)} {run.records === 1 ? 'record' : 'records'} ·{' '}
              {fmt(run.bases)} bases
              {preview.cut && ' · showing the start; copy or download for all'}
            </span>
            <button type="button" className="rg-btn rg-btn--sm" onClick={download}>
              <Icon name="download" /> Download .fa
            </button>
          </div>
          <Copyable
            variant="block"
            value={run.fasta}
            label="Copy the FASTA"
            className="rg-region__fasta"
          >
            {preview.text}
          </Copyable>
        </div>
      )}
    </section>
  );
}

/** The start of a FASTA, cut at a line end. */
function previewOf(fasta: string): { text: string; cut: boolean } {
  if (fasta.length <= PREVIEW_CHARS) return { text: fasta, cut: false };
  const end = fasta.lastIndexOf('\n', PREVIEW_CHARS);
  return { text: `${fasta.slice(0, end > 0 ? end : PREVIEW_CHARS)}\n…`, cut: true };
}

function shorten(line: string, max = 60): string {
  return line.length > max ? `${line.slice(0, max)}…` : line;
}
