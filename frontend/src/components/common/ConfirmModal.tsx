/**
 * The one confirmation dialog. Every destructive action in the app goes
 * through it; no component writes modal JSX of its own (that is the specific
 * mistake the refget frontend made three times).
 */

import { useEffect, useState } from 'react';
import type { ReactNode } from 'react';
import { BaseModal } from './BaseModal';
import { cn } from '../../utils/cn';
import { ApiError } from '../../services/http';

export interface ConfirmModalProps {
  isOpen: boolean;
  onClose: () => void;
  title: string;
  body: ReactNode;
  confirmLabel: string;
  destructive?: boolean;
  /**
   * When set, the confirm button stays disabled until the user types this
   * exactly. Reserved for cascades — deleting a genome takes its aliases,
   * asset groups and assets with it, unconditionally.
   */
  confirmationText?: string;
  onConfirm: () => Promise<unknown> | unknown;
}

export function ConfirmModal({
  isOpen,
  onClose,
  title,
  body,
  confirmLabel,
  destructive = false,
  confirmationText,
  onConfirm,
}: ConfirmModalProps) {
  const [typed, setTyped] = useState('');
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (isOpen) {
      setTyped('');
      setBusy(false);
      setError(null);
    }
  }, [isOpen]);

  const gated = confirmationText !== undefined && typed !== confirmationText;

  const handleConfirm = async () => {
    setBusy(true);
    setError(null);
    try {
      await onConfirm();
      setBusy(false);
      onClose();
    } catch (caught) {
      setBusy(false);
      setError(caught instanceof ApiError ? caught.detail : String(caught));
    }
  };

  return (
    <BaseModal isOpen={isOpen} onClose={onClose} title={title} size="md">
      <div className="flex flex-col gap-4">
        <div className="text-sm">{body}</div>

        {confirmationText !== undefined && (
          <div className="rg-field">
            <label className="rg-field__label" htmlFor="confirm-modal-phrase">
              Type <code className="rg-code rg-code--inline">{confirmationText}</code> to
              confirm
            </label>
            <input
              id="confirm-modal-phrase"
              className="rg-field__input rg-field__input--mono"
              type="text"
              value={typed}
              autoComplete="off"
              onChange={(event) => setTyped(event.target.value)}
            />
          </div>
        )}

        {error && (
          <p className="rg-field__error" role="alert">
            {error}
          </p>
        )}

        <BaseModal.Footer>
          <button type="button" className="rg-btn" onClick={onClose} disabled={busy}>
            Cancel
          </button>
          <button
            type="button"
            className={cn('rg-btn', destructive ? 'rg-btn--danger' : 'rg-btn--primary')}
            onClick={() => void handleConfirm()}
            disabled={busy || gated}
          >
            {busy ? 'Working…' : confirmLabel}
          </button>
        </BaseModal.Footer>
      </div>
    </BaseModal>
  );
}
