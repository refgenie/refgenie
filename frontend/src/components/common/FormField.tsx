import { Children, cloneElement, isValidElement } from 'react';
import type { ReactNode } from 'react';
import { cn } from '../../utils/cn';
import { useHelpDisclosure } from './useHelpDisclosure';

export interface FormFieldProps {
  /** Must match the control's `id`, so the label actually labels it. */
  htmlFor: string;
  label: string;
  /** Help text, shown behind a "?" toggle and always linked via `aria-describedby`. */
  hint?: ReactNode;
  /** Field error, shown inline (never as a toast) and never collapsed. */
  error?: string;
  required?: boolean;
  children: ReactNode;
}

type DescribableProps = { id?: string; 'aria-describedby'?: string };

/** Add the hint id to `aria-describedby` on the child whose `id` matches `htmlFor`. */
function describeControl(children: ReactNode, controlId: string, hintId: string): ReactNode {
  return Children.map(children, (child) => {
    if (!isValidElement<DescribableProps>(child)) return child;
    if (child.props.id !== controlId) return child;
    const existing = child.props['aria-describedby'];
    return cloneElement(child, {
      'aria-describedby': existing ? `${existing} ${hintId}` : hintId,
    });
  });
}

export function FormField({
  htmlFor,
  label,
  hint,
  error,
  required,
  children,
}: FormFieldProps) {
  const hintId = `${htmlFor}-hint`;
  const help = useHelpDisclosure(hintId, label);

  return (
    <div className={cn('rg-field', error && 'rg-field--invalid')}>
      <div className="rg-field__label-cell">
        <label className="rg-field__label" htmlFor={htmlFor}>
          {label}
          {required && (
            <span className="rg-field__required" aria-hidden="true">
              {' '}
              *
            </span>
          )}
        </label>
        {hint && <button {...help.toggleProps}>?</button>}
      </div>

      <div className="rg-field__control">
        {hint ? describeControl(children, htmlFor, hintId) : children}
      </div>

      {hint && <p {...help.hintProps}>{hint}</p>}

      {error && (
        <p className="rg-field__error" role="alert">
          {error}
        </p>
      )}
    </div>
  );
}
