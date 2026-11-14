import { useCallback, useEffect, useId, useRef, useState } from 'react';
import { useSearchParams } from 'react-router-dom';
import { FastaResult } from './FastaResult';
import { extractRegion } from '../../services/regionExtract';
import {
  MAX_EXTRACT_BASES,
  exampleRegion,
  parseRegion,
  regionFasta,
  regionFileName,
  type ParsedRegion,
} from '../../utils/parseRegion';
import type { SequenceResponse } from '../../types/api';

const REGION_PARAM = 'region';

export interface RegionExtractProps {
  storeUrl: string;
  /** All of the genome's sequences, to match names and aliases against. */
  sequences: SequenceResponse[];
  value: string;
  onChange: (value: string) => void;
}

type Result =
  | { kind: 'idle' }
  | { kind: 'working' }
  | { kind: 'error'; message: string }
  | { kind: 'done'; region: ParsedRegion; fasta: string; length: number };

/**
 * The Extract tab's region box: type `chr1:1,000,000-1,000,500`, get the
 * bases as FASTA. Only the bytes for the region are downloaded. The region is
 * kept in the URL (`?region=`), so a link to it opens with the bases shown.
 */
export function RegionExtract({
  storeUrl,
  sequences,
  value,
  onChange,
}: RegionExtractProps) {
  const inputId = useId();
  const hintId = useId();
  const [params, setParams] = useSearchParams();
  const [result, setResult] = useState<Result>({ kind: 'idle' });
  const controller = useRef<AbortController | null>(null);
  // The last region run, input or normalized, so the URL it writes does not
  // run it again.
  const lastRun = useRef<string | null>(null);

  const run = useCallback(
    async (text: string) => {
      controller.current?.abort();
      lastRun.current = text;
      let region: ParsedRegion;
      try {
        region = parseRegion(text, sequences);
      } catch (error) {
        setResult({ kind: 'error', message: messageOf(error) });
        return;
      }
      const seq = sequences.find((row) => row.name === region.name) as SequenceResponse;
      lastRun.current = region.label;
      setParams(
        (previous) => {
          const next = new URLSearchParams(previous);
          next.set(REGION_PARAM, region.label);
          return next;
        },
        { replace: true },
      );
      const abort = new AbortController();
      controller.current = abort;
      setResult({ kind: 'working' });
      try {
        const bases = await extractRegion({
          baseUrl: storeUrl,
          seq,
          start: region.start,
          end: region.end,
          signal: abort.signal,
        });
        if (abort.signal.aborted) return;
        setResult({
          kind: 'done',
          region,
          fasta: regionFasta(region.label, bases),
          length: bases.length,
        });
      } catch (error) {
        if (abort.signal.aborted) return;
        setResult({ kind: 'error', message: messageOf(error) });
      }
    },
    [sequences, storeUrl, setParams],
  );

  // A `?region=` link (on arrival, or back/forward to another one) fills the
  // box and runs.
  const linked = params.get(REGION_PARAM);
  useEffect(() => {
    if (!linked || linked === lastRun.current) return;
    onChange(linked);
    void run(linked);
  }, [linked, onChange, run]);

  // Leaving stops a read in flight; forgetting it lets a remount run it again.
  useEffect(
    () => () => {
      controller.current?.abort();
      lastRun.current = null;
    },
    [],
  );

  const example = exampleRegion(sequences);

  return (
    <section className="rg-region" aria-label="Extract a region">
      <form
        className="rg-region__form"
        onSubmit={(event) => {
          event.preventDefault();
          void run(value);
        }}
      >
        <label className="rg-region__label" htmlFor={inputId}>
          Region
        </label>
        <input
          id={inputId}
          className="rg-search__input rg-region__input"
          type="text"
          spellCheck={false}
          autoComplete="off"
          placeholder={example ?? 'name:start-end'}
          aria-describedby={hintId}
          value={value}
          onChange={(event) => onChange(event.target.value)}
        />
        <button type="submit" className="rg-btn" disabled={result.kind === 'working'}>
          {result.kind === 'working' ? 'Extracting…' : 'Extract'}
        </button>
        {example && (
          <button
            type="button"
            className="rg-btn"
            disabled={result.kind === 'working'}
            onClick={() => {
              onChange(example);
              void run(example);
            }}
          >
            Try an example
          </button>
        )}
      </form>
      <p id={hintId} className="rg-region__hint">
        1-based, inclusive, like samtools faidx. Up to {MAX_EXTRACT_BASES.toLocaleString()}{' '}
        bases; only those bytes are downloaded.
      </p>

      {result.kind === 'error' && (
        <p className="rg-region__error" role="alert">
          {result.message}
        </p>
      )}
      {result.kind === 'done' && (
        <FastaResult
          fasta={result.fasta}
          fileName={regionFileName(result.region.label)}
          meta={
            `${result.region.label} · ${result.length.toLocaleString()} ` +
            (result.length === 1 ? 'base' : 'bases') +
            (result.region.clamped ? ' (cut at the end of the sequence)' : '')
          }
        />
      )}
    </section>
  );
}

function messageOf(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}
