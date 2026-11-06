import { useContext } from 'react';
import { ToastContext } from '../components/common/toastContext';
import type { ToastApi } from '../components/common/toastContext';

/** Toasts for short mutations and "queued" notices; field errors render inline instead. */
export function useToast(): ToastApi {
  const api = useContext(ToastContext);
  if (!api) throw new Error('useToast must be used inside <ToastProvider>');
  return api;
}
