import type { ReactNode } from 'react';
import { Icon } from '../common/Icon';
import { Copyable } from '../common/Copyable';
import { downloadBlob } from '../../services/saveStream';

/** How much of the FASTA to show on the page; copy and download get it all. */
const PREVIEW_CHARS = 20_000;

/** An extracted FASTA: a line about it, a download button, and the text to copy. */
export function FastaResult({
  fasta,
  fileName,
  meta,
}: {
  fasta: string;
  fileName: string;
  meta: ReactNode;
}) {
  const preview = previewOf(fasta);
  return (
    <div className="rg-region__out">
      <div className="rg-region__bar">
        <span className="rg-region__meta">
          {meta}
          {preview.cut && ' · showing the start; copy or download for all'}
        </span>
        <button
          type="button"
          className="rg-btn rg-btn--sm"
          onClick={() => downloadBlob(fileName, new Blob([fasta], { type: 'text/plain' }))}
        >
          <Icon name="download" /> Download .fa
        </button>
      </div>
      <Copyable
        variant="block"
        value={fasta}
        label="Copy the FASTA"
        className="rg-region__fasta"
      >
        {preview.text}
      </Copyable>
    </div>
  );
}

/** The start of a FASTA, cut at a line end. */
function previewOf(fasta: string): { text: string; cut: boolean } {
  if (fasta.length <= PREVIEW_CHARS) return { text: fasta, cut: false };
  const end = fasta.lastIndexOf('\n', PREVIEW_CHARS);
  return { text: `${fasta.slice(0, end > 0 ? end : PREVIEW_CHARS)}\n…`, cut: true };
}
