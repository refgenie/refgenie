import { useState } from 'react';
import { BaseModal } from './BaseModal';

export interface JsonBlockProps {
  value: unknown;
  /** The button text. */
  label?: string;
  /** The modal title; defaults to `label`. */
  title?: string;
}

/** A button that opens raw JSON in a modal. No dangerouslySetInnerHTML anywhere. */
export function JsonBlock({ value, label = 'Raw JSON', title }: JsonBlockProps) {
  const [open, setOpen] = useState(false);
  const heading = title ?? label;
  return (
    <>
      <button
        type="button"
        className="rg-btn rg-btn--sm rg-btn--bare"
        aria-haspopup="dialog"
        onClick={() => setOpen(true)}
      >
        {label}
      </button>
      <BaseModal isOpen={open} onClose={() => setOpen(false)} title={heading} size="lg">
        <pre className="rg-code" aria-label={heading}>
          {JSON.stringify(value, null, 2)}
        </pre>
      </BaseModal>
    </>
  );
}
