import type { ReactNode } from 'react';
import { cn } from '../../utils/cn';

export interface FormFieldProps {
  /** Must match the control's `id`, so the label actually labels it. */
  htmlFor: string;
  label: string;
  hint?: ReactNode;
  /**
   * Field-level errors render HERE, never as a toast: a toast covering the
   * field it is about is worse than no message.
   */
  error?: string;
  required?: boolean;
  children: ReactNode;
}

export function FormField({
  htmlFor,
  label,
  hint,
  error,
  required,
  children,
}: FormFieldProps) {
  return (
    <div className={cn('rg-field', error && 'rg-field--invalid')}>
      <label className="rg-field__label" htmlFor={htmlFor}>
        {label}
        {required && (
          <span className="rg-field__required" aria-hidden="true">
            {' '}
            *
          </span>
        )}
      </label>
      {children}
      {hint && <p className="rg-field__hint">{hint}</p>}
      {error && (
        <p className="rg-field__error" role="alert">
          {error}
        </p>
      )}
    </div>
  );
}
