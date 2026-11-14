import type { ReactNode } from 'react';
import { CopyButton } from './CopyButton';
import { cn } from '../../utils/cn';

export interface CopyableProps {
  /** The exact text the button copies. */
  value: string;
  /**
   * `inline`: a short value (an alias, a digest) with its button right after it.
   * `block`: a code snippet, with its button in the block's top-right corner.
   */
  variant?: 'inline' | 'block';
  /** The button's accessible name. */
  label?: string;
  /** What to show, when it differs from `value` (a shortened digest, say). */
  children?: ReactNode;
  /** Monospace for the inline variant. Blocks are always monospace. */
  mono?: boolean;
  className?: string;
}

/** A value and its copy button, side by side. Every copyable thing renders through this. */
export function Copyable({
  value,
  variant = 'inline',
  label,
  children,
  mono = false,
  className,
}: CopyableProps) {
  const button = (
    <CopyButton
      value={value}
      label={label ?? (variant === 'block' ? 'Copy' : `Copy ${value}`)}
    />
  );
  if (variant === 'block') {
    return (
      <div className={cn('rg-copyable-block', className)}>
        <pre className="rg-code rg-copyable-block__code">
          {children ?? value}
        </pre>
        {button}
      </div>
    );
  }
  return (
    <span className={cn('rg-copyable', className)}>
      <span className={cn('rg-copyable__text', mono && 'rg-copyable__text--mono')}>
        {children ?? value}
      </span>
      {button}
    </span>
  );
}

/** A shell command to copy, shown as inline code. */
export function Command({ value }: { value: string }) {
  return (
    <Copyable value={value} label={`Copy: ${value}`} className="rg-command">
      <code className="rg-code rg-code--inline">{value}</code>
    </Copyable>
  );
}
