import { Children, cloneElement, isValidElement } from 'react';
import type { ReactNode } from 'react';
import { cn } from '../../utils/cn';
import { useHelpDisclosure } from './useHelpDisclosure';

export interface FormFieldProps {
  /** Must match the control's `id`, so the label actually labels it. */
  htmlFor: string;
  label: string;
  /**
   * The explainer. It renders behind a "?" toggle rather than under the
   * control: a permanent third line per field is what made these forms scroll.
   * Collapsed is a *visual* state only — the text stays in the DOM and stays
   * wired to the control with `aria-describedby`, so a screen reader reads it
   * either way.
   */
  hint?: ReactNode;
  /**
   * Field-level errors render HERE, never as a toast: a toast covering the
   * field it is about is worse than no message. They are also never collapsed
   * behind the help toggle — an error nobody can see is not a message.
   */
  error?: string;
  required?: boolean;
  children: ReactNode;
}

type DescribableProps = { id?: string; 'aria-describedby'?: string };

/**
 * Point the control at the hint without asking eight call sites to repeat the
 * id. The control is the child whose `id` matches `htmlFor` — the same
 * contract the `<label for>` already depends on — so a sibling that is not the
 * control (GenomeSelect's `<datalist>`) is left alone.
 */
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
