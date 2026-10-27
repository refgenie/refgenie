import { useContext } from 'react';
import { ToastContext } from '../components/common/toastContext';
import type { ToastApi } from '../components/common/toastContext';

/**
 * Toasts are for short synchronous mutations and for the one-line "queued"
 * acknowledgement of a long operation. Never for a field-level validation
 * error: a toast hides the field it is about.
 */
export function useToast(): ToastApi {
  const api = useContext(ToastContext);
  if (!api) throw new Error('useToast must be used inside <ToastProvider>');
  return api;
}
