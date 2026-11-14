import { useId } from 'react';
import { CodeSnippets } from '../common/CodeSnippets';
import { SequencePicker } from './SequencePicker';
import { downloadSnippets } from '../../utils/codeSnippets';
import { installCommand } from '../../utils/install';
import { formatBytes } from '../../utils/format';
import { useUiConfig } from '../../hooks/useUiConfig';
import type { FastaDownloadState } from '../../hooks/useFastaDownload';

export interface GenomeDownloadPanelProps {
  digest: string;
  /** The genome's refget store. */
  storeUrl: string;
  /** A name for the file and for refgenie commands; the digest if none. */
  alias: string | undefined;
  /** The in-browser download. */
  download: FastaDownloadState;
  /** True if the genome has a built fasta asset for refgenie to pull. */
  hasFastaAsset: boolean;
}

/**
 * The genome page's Download tab: the FASTA rebuilt in the browser from the
 * genome's refget store, or with a command-line tool. The download state
 * lives on the page, so switching tabs does not stop a download.
 */
export function GenomeDownloadPanel({
  digest,
  storeUrl,
  alias,
  download,
  hasFastaAsset,
}: GenomeDownloadPanelProps) {
  const version = useUiConfig().refgenie_version;
  const choiceName = useId();
  const { collection, status, choosing, selected, count, totalBytes, estimatedBytes } =
    download;
  const saving = status.kind === 'saving';
  return (
    <div>
      <section className="rg-flat-section">
        <h2>In your browser</h2>
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
              `about ${formatBytes(estimatedBytes)} compressed ` +
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
          </button>
          {saving && (
            <>
              <span role="status">
                Saving… {Math.floor(status.fraction * 100)}% · {formatBytes(status.written)}
              </span>
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
      </section>
      <section className="rg-flat-section">
        <h2>With other tools</h2>
        <p className="rg-muted text-sm mb-3">
          The same FASTA, downloaded on your own computer or server.
        </p>
        <CodeSnippets
          idPrefix="download-code"
          label="Tools"
          snippets={downloadSnippets({
            digest,
            storeUrl,
            alias,
            names: download.names,
            hasFastaAsset,
            refgenieInstall: installCommand(version),
          })}
        />
      </section>
    </div>
  );
}
