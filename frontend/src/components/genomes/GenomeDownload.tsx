import { CodeSnippets } from '../common/CodeSnippets';
import { FastaDownloadPanel } from './FastaDownload';
import { downloadSnippets } from '../../utils/codeSnippets';
import type { FastaDownloadState } from '../../hooks/useFastaDownload';

export interface GenomeDownloadPanelProps {
  digest: string;
  /** The genome's refget store. */
  storeUrl: string;
  /** A name for the file and for refgenie commands; the digest if none. */
  alias: string | undefined;
  /** The in-browser download; null when the browser may not read the store. */
  download: FastaDownloadState | null;
  /** True if the genome has a built fasta asset for refgenie to pull. */
  hasFastaAsset: boolean;
}

/** The genome page's Download tab: the FASTA in the browser, or with a command-line tool. */
export function GenomeDownloadPanel({
  digest,
  storeUrl,
  alias,
  download,
  hasFastaAsset,
}: GenomeDownloadPanelProps) {
  return (
    <div>
      {download && (
        <section className="rg-flat-section">
          <h2>In your browser</h2>
          <FastaDownloadPanel download={download} />
        </section>
      )}
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
            names: download?.names,
            hasFastaAsset,
          })}
        />
      </section>
    </div>
  );
}
