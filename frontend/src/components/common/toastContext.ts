/**
 * Toast context, split from the provider component so the module exports only
 * values (the provider file stays components-only for react-refresh).
 */

import { createContext } from 'react';

export type ToastTone = 'success' | 'error' | 'info';

export interface ToastAction {
  label: string;
  onClick: () => void;
}

export interface ToastOptions {
  tone?: ToastTone;
  /** Overrides the tone default: 5s for success/info, sticky for error. */
  durationMs?: number | null;
  action?: ToastAction;
}

export interface ToastRecord extends Required<Pick<ToastOptions, 'tone'>> {
  id: number;
  message: string;
  action?: ToastAction;
}

export interface ToastApi {
  show: (message: string, options?: ToastOptions) => number;
  success: (message: string, options?: Omit<ToastOptions, 'tone'>) => number;
  error: (message: string, options?: Omit<ToastOptions, 'tone'>) => number;
  info: (message: string, options?: Omit<ToastOptions, 'tone'>) => number;
  dismiss: (id: number) => void;
}

/** Errors are sticky: a failure the user did not read is a failure lost. */
export const TOAST_DEFAULT_DURATION_MS = 5000;

export const ToastContext = createContext<ToastApi | undefined>(undefined);
