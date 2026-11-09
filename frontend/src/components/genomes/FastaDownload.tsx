import { useId } from 'react';
import { SequencePicker } from './SequencePicker';
import { formatBytes } from '../../utils/format';
import { approxBytes, savingText } from '../../utils/fastaDownload';
import type { FastaDownloadState } from '../../hooks/useFastaDownload';

export interface FastaDownloadPanelProps {
  download: FastaDownloadState;
}

/**
 * The genome's in-browser FASTA download, on the Download tab. The FASTA is
 * rebuilt in the browser from the genome's refget store, so it belongs to the
 * genome itself, not to any asset. The download state lives on the page, so
 * switching tabs does not stop a download.
 */
export function FastaDownloadPanel({ download }: FastaDownloadPanelProps) {
  const choiceName = useId();
  const { collection, status, choosing, selected, count, totalBytes, estimatedBytes } =
    download;
  const saving = status.kind === 'saving';

  return (
    <div>
      <fieldset className="rg-fasta-choice" disabled={saving}>
        <legend className="sr-only">Sequences to include</legend>
        <label className="flex items-center gap-2">
          <input
            type="radio"
            name={choiceName}
            checked={!choosing}
            onChange={() => download.setChoosing(false)}
          />
          All sequences
        </label>
        <label className="flex items-center gap-2">
          <input
            type="radio"
            name={choiceName}
            checked={choosing}
            disabled={!collection}
            onChange={() => download.setChoosing(true)}
          />
          Choose sequences…
        </label>
      </fieldset>
      {choosing && collection && (
        <SequencePicker
          sequences={collection.seqs}
          selected={selected}
          onChange={download.setSelected}
          unavailable={download.unavailable}
        />
      )}
      <p className="rg-muted text-sm mt-4">
        {count === null || totalBytes === null || estimatedBytes === null
          ? 'Working out the size…'
          : `${count.toLocaleString()} ${count === 1 ? 'sequence' : 'sequences'}, ` +
            `about ${approxBytes(estimatedBytes)} compressed ` +
            `(${formatBytes(totalBytes)} uncompressed). `}
        Saved as a bgzip-compressed .fa.gz, built in your browser from the refget store.
      </p>
      <div className="flex flex-wrap items-center gap-4 mt-4">
        <button
          type="button"
          className="rg-btn rg-btn--primary"
          onClick={download.start}
          disabled={!collection || count === 0 || saving}
        >
          Download .fa.gz
          {/* Fixed width, so the size arriving moves nothing. */}
          <span className="rg-size-slot">
            {estimatedBytes !== null ? (
              count === 0 ? null : ` · about ${approxBytes(estimatedBytes)}`
            ) : (
              <span className="rg-skeleton rg-skeleton--inline rg-size-slot__ph" />
            )}
          </span>
        </button>
        {saving && (
          <>
            <span role="status">{savingText(status.fraction, status.written)}</span>
            <button type="button" className="rg-btn rg-btn--sm" onClick={download.cancel}>
              Cancel
            </button>
          </>
        )}
        {status.kind === 'saved' && <span role="status">Saved</span>}
      </div>
      {status.kind === 'error' && (
        <div className="rg-state rg-state--error mt-4" role="alert">
          {status.message}
        </div>
      )}
    </div>
  );
}
