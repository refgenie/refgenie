import { useEffect, useMemo, useRef, useState } from 'react';
import { useSearchParams } from 'react-router-dom';
import {
  useGenomeComparison,
  useGenomeSequences,
} from '../../hooks/queries/useGenomeSequences';
import type { CompareCandidate } from '../../utils/compareCandidates';
import { ErrorState, LoadingState } from '../common/states';
import { CodeSnippets } from '../common/CodeSnippets';
import { CompareFasta } from './CompareFasta';
import { CompareGenomes } from './CompareGenomes';
import { SequenceStats } from './SequenceStats';
import { SequenceTable } from './SequenceTable';
import { RegionExtract } from './RegionExtract';
import { BedExtract } from './BedExtract';
import { defaultRegion, parseRegion, type ParsedRegion } from '../../utils/parseRegion';
import { bedSnippets, regionSnippets } from '../../utils/codeSnippets';
import { cn } from '../../utils/cn';
import { extractHref } from '../../utils/genomePage';
import {
  isExtractable,
  sequenceFileUrl,
  sequenceTemplate,
} from '../../services/regionExtract';
import {
  DEFAULT_COLLECTION_TEMPLATE,
  DEFAULT_SEQ_TEMPLATE,
  seqUrl,
} from '../../services/fastaExport';
import type { SequenceRow } from './types';

export interface GenomePanelProps {
  digest: string;
}

/**
 * The genome page's Sequences tab: a one-line summary and the sequence table.
 * `canExtract` adds a per-row link to the Extract tab.
 */
export function SequencesPanel({
  digest,
  canExtract = false,
}: GenomePanelProps & { canExtract?: boolean }) {
  const sequences = useGenomeSequences(digest);
  const data = sequences.data;

  if (sequences.error) {
    return <ErrorState error={sequences.error} onRetry={() => sequences.refetch()} />;
  }
  return (
    <div className="flex flex-col gap-4">
      <h2 className="sr-only">Sequences</h2>
      <p className="rg-muted text-sm">
        Digests (SHA-512/24u, the GA4GH refget ID, and MD5) are IDs computed from a
        sequence&apos;s bases: the same digest means exactly the same sequence. Alphabet:
        dna2bit holds only A, C, G and T; dna3bit adds N; dnaio allows any IUPAC letter.
      </p>
      {data ? (
        <SequenceTable
          sequences={data.sequences}
          extractHref={canExtract ? extractHref : undefined}
          aside={<SequenceStats sequences={data.sequences} />}
        />
      ) : (
        <LoadingState rows={12} label="Loading sequences" />
      )}
    </div>
  );
}

const NO_LIMIT = Number.POSITIVE_INFINITY;

function tryParse(text: string, sequences: SequenceRow[]): ParsedRegion | null {
  try {
    return parseRegion(text, sequences, NO_LIMIT);
  } catch {
    return null;
  }
}

type ExtractMode = 'one' | 'bed';

const MODES: { id: ExtractMode; label: string }[] = [
  { id: 'one', label: 'One region' },
  { id: 'bed', label: 'BED file' },
];

/**
 * The genome page's Extract tab: one region, or every region of a BED file.
 * Each reads the bases in the browser (when `browserStoreUrl` is set), and
 * shows the same in code, for each tool that can read the genome's refget
 * store.
 */
export function ExtractPanel({
  digest,
  storeUrl,
  browserStoreUrl,
}: GenomePanelProps & {
  /** The genome's refget store, named in the code snippets. */
  storeUrl: string;
  /** Set when the browser may read the store itself. */
  browserStoreUrl?: string | null;
}) {
  const sequences = useGenomeSequences(digest);
  const data = sequences.data;
  const [params] = useSearchParams();
  const [region, setRegion] = useState(() => params.get('region') ?? '');
  const regionInput = useRef<HTMLInputElement | null>(null);
  const [mode, setMode] = useState<ExtractMode>('one');

  // The snippets follow the box: the region typed, else the last one that
  // parsed, else the first readable sequence's first 1,000 bases.
  const parsed = useMemo(
    () => (data ? tryParse(region, data.sequences) : null),
    [region, data],
  );
  const [lastGood, setLastGood] = useState<ParsedRegion | null>(null);
  if (parsed && parsed.label !== lastGood?.label) setLastGood(parsed);
  const fallback = useMemo(() => {
    const first = data?.sequences.find((row) => isExtractable(row.alphabet));
    return first && data ? tryParse(defaultRegion(first), data.sequences) : null;
  }, [data]);
  const shown = parsed ?? lastGood ?? fallback;
  const seq = shown ? data?.sequences.find((row) => row.name === shown.name) : undefined;

  // The `.seq` URL by the store's own path template, when the browser may read
  // the store's manifest; else the usual template.
  const seqDigest = seq?.sha512t24u;
  const [fileUrls, setFileUrls] = useState<Record<string, string>>({});
  useEffect(() => {
    if (!browserStoreUrl || !seqDigest) return;
    let stale = false;
    sequenceFileUrl(browserStoreUrl, seqDigest)
      .then((url) => {
        if (!stale) setFileUrls((known) => ({ ...known, [seqDigest]: url }));
      })
      .catch(() => {});
    return () => {
      stale = true;
    };
  }, [browserStoreUrl, seqDigest]);

  // The store's own `.seq` template, for the BED code; else the usual one.
  const [seqTemplate, setSeqTemplate] = useState(DEFAULT_SEQ_TEMPLATE);
  useEffect(() => {
    if (!browserStoreUrl || mode !== 'bed') return;
    let stale = false;
    sequenceTemplate(browserStoreUrl)
      .then((template) => {
        if (!stale) setSeqTemplate(template);
      })
      .catch(() => {});
    return () => {
      stale = true;
    };
  }, [browserStoreUrl, mode]);

  if (sequences.error) {
    return <ErrorState error={sequences.error} onRetry={() => sequences.refetch()} />;
  }
  if (!data) return <LoadingState rows={8} label="Loading sequences" />;
  return (
    <div>
      <div className="rg-segmented" role="radiogroup" aria-label="What to extract">
        {MODES.map(({ id, label }) => (
          <button
            key={id}
            type="button"
            role="radio"
            aria-checked={mode === id}
            className={cn('rg-segmented__option', mode === id && 'rg-segmented__option--on')}
            onClick={() => setMode(id)}
          >
            {label}
          </button>
        ))}
      </div>
      {browserStoreUrl && (
        <section className="rg-flat-section">
          <h2>In your browser</h2>
          {/* Both stay mounted, so switching keeps each one's input and result. */}
          <div hidden={mode !== 'one'}>
            <RegionExtract
              storeUrl={browserStoreUrl}
              sequences={data.sequences}
              value={region}
              onChange={setRegion}
              inputRef={regionInput}
            />
          </div>
          <div hidden={mode !== 'bed'}>
            <BedExtract storeUrl={browserStoreUrl} sequences={data.sequences} />
          </div>
        </section>
      )}
      {mode === 'bed' && (
        <section className="rg-flat-section">
          <h2>In your own code</h2>
          <p className="rg-muted text-sm mb-3">
            Every region of a BED file (saved as{' '}
            <code className="rg-code rg-code--inline">regions.bed</code>) from this
            genome&apos;s refget store, into a FASTA file, with your own tools. No size
            limit.
          </p>
          <CodeSnippets
            idPrefix="extract-bed-code"
            label="Tools"
            snippets={bedSnippets({
              digest,
              storeUrl,
              seqTemplate,
              collectionTemplate: DEFAULT_COLLECTION_TEMPLATE,
            })}
          />
        </section>
      )}
      {mode === 'one' && shown && seq && (
        <section className="rg-flat-section">
          <h2>In your own code</h2>
          <p className="rg-muted text-sm mb-3">
            The region <code className="rg-code rg-code--inline">{shown.label}</code>{' '}
            from this genome&apos;s refget store, with your own tools.
          </p>
          <CodeSnippets
            idPrefix="extract-code"
            label="Tools"
            snippets={regionSnippets({
              digest,
              storeUrl,
              region: shown,
              seq: { sha512t24u: seq.sha512t24u, alphabet: seq.alphabet ?? '' },
              seqFileUrl:
                fileUrls[seq.sha512t24u] ??
                seqUrl(storeUrl, DEFAULT_SEQ_TEMPLATE, seq.sha512t24u),
            })}
          />
        </section>
      )}
    </div>
  );
}

/**
 * The genome page's Compare tab: compare with another genome of the species,
 * or (with `canCompareFile`) with the user's own FASTA file.
 */
export function ComparePanel({
  digest,
  candidates,
  canCompareFile = false,
}: GenomePanelProps & { candidates: CompareCandidate[] | null; canCompareFile?: boolean }) {
  const sequences = useGenomeSequences(digest);
  const [compareWith, setCompareWith] = useState<string>('');
  const comparison = useGenomeComparison(digest, compareWith || undefined);

  if (sequences.error) {
    return <ErrorState error={sequences.error} onRetry={() => sequences.refetch()} />;
  }
  const data = sequences.data;
  return (
    <div>
      <CompareGenomes
        digest={digest}
        candidates={candidates ?? []}
        comparison={comparison.data}
        loading={comparison.isPending}
        error={comparison.error}
        onSelect={setCompareWith}
        selected={compareWith || undefined}
      />
      {canCompareFile && <CompareFasta digest={digest} sequences={data?.sequences} />}
    </div>
  );
}
