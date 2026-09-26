/**
 * In-house toasts, ~80 lines, rather than `react-hot-toast`: that library
 * injects inline styles, which the locked styling decision forbids.
 *
 * Routing rules this implements (see the feedback conventions):
 *  - field-level validation errors never come here, they render inline;
 *  - short synchronous mutations toast and refetch;
 *  - long operations toast once ("queued") and then live in the job console.
 */

import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import type { ReactNode } from 'react';
import { cn } from '../../utils/cn';
import {
  TOAST_DEFAULT_DURATION_MS,
  ToastContext,
} from './toastContext';
import type { ToastApi, ToastOptions, ToastRecord } from './toastContext';

export function ToastProvider({ children }: { children: ReactNode }) {
  const [toasts, setToasts] = useState<ToastRecord[]>([]);
  const nextId = useRef(1);
  const timers = useRef(new Map<number, ReturnType<typeof setTimeout>>());

  const dismiss = useCallback((id: number) => {
    const timer = timers.current.get(id);
    if (timer !== undefined) {
      clearTimeout(timer);
      timers.current.delete(id);
    }
    setToasts((current) => current.filter((toast) => toast.id !== id));
  }, []);

  const show = useCallback(
    (message: string, options: ToastOptions = {}) => {
      const tone = options.tone ?? 'info';
      const id = nextId.current;
      nextId.current += 1;
      setToasts((current) => [...current, { id, message, tone, action: options.action }]);

      const duration =
        options.durationMs === undefined
          ? tone === 'error'
            ? null
            : TOAST_DEFAULT_DURATION_MS
          : options.durationMs;
      if (duration !== null) {
        timers.current.set(
          id,
          setTimeout(() => dismiss(id), duration),
        );
      }
      return id;
    },
    [dismiss],
  );

  useEffect(() => {
    const pending = timers.current;
    return () => {
      for (const timer of pending.values()) clearTimeout(timer);
      pending.clear();
    };
  }, []);

  const api = useMemo<ToastApi>(
    () => ({
      show,
      dismiss,
      success: (message, options) => show(message, { ...options, tone: 'success' }),
      error: (message, options) => show(message, { ...options, tone: 'error' }),
      info: (message, options) => show(message, { ...options, tone: 'info' }),
    }),
    [show, dismiss],
  );

  return (
    <ToastContext.Provider value={api}>
      {children}
      <div className="rg-toast-region" role="status" aria-live="polite">
        {toasts.map((toast) => (
          <div className={cn('rg-toast', `rg-toast--${toast.tone}`)} key={toast.id}>
            <span className="rg-toast__message">{toast.message}</span>
            {toast.action && (
              <button
                type="button"
                className="rg-btn rg-btn--sm"
                onClick={() => {
                  toast.action?.onClick();
                  dismiss(toast.id);
                }}
              >
                {toast.action.label}
              </button>
            )}
            <button
              type="button"
              className="rg-toast__close"
              aria-label="Dismiss notification"
              onClick={() => dismiss(toast.id)}
            >
              &times;
            </button>
          </div>
        ))}
      </div>
    </ToastContext.Provider>
  );
}
