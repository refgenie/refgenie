/** The shared confirmation dialog for destructive actions. */

import { useEffect, useState } from 'react';
import type { ReactNode } from 'react';
import { BaseModal } from './BaseModal';
import { cn } from '../../utils/cn';
import { useAction } from '../../hooks/useAction';

export interface ConfirmModalProps {
  isOpen: boolean;
  onClose: () => void;
  title: string;
  body: ReactNode;
  confirmLabel: string;
  destructive?: boolean;
  /** If set, confirm stays disabled until the user types this. Use for cascading deletes. */
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
  const { busy, error, setError, run } = useAction();

  useEffect(() => {
    if (isOpen) {
      setTyped('');
      setError(null);
    }
  }, [isOpen, setError]);

  const gated = confirmationText !== undefined && typed !== confirmationText;

  const handleConfirm = () =>
    run(async () => {
      await onConfirm();
      onClose();
    });

  return (
    <BaseModal isOpen={isOpen} onClose={onClose} title={title} size="md">
      <div className="flex flex-col gap-4">
        <div className="text-sm">{body}</div>

        {confirmationText !== undefined && (
          <div className="rg-field rg-field--stacked">
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
