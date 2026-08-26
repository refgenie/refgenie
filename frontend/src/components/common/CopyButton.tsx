import { useCallback, useEffect, useRef, useState } from 'react';

export interface CopyButtonProps {
  value: string;
  label?: string;
}

export function CopyButton({ value, label = 'Copy' }: CopyButtonProps) {
  const [copied, setCopied] = useState(false);
  const timer = useRef<ReturnType<typeof setTimeout> | undefined>(undefined);

  useEffect(() => () => clearTimeout(timer.current), []);

  const copy = useCallback(async () => {
    try {
      await navigator.clipboard.writeText(value);
      setCopied(true);
      clearTimeout(timer.current);
      timer.current = setTimeout(() => setCopied(false), 2000);
    } catch {
      setCopied(false);
    }
  }, [value]);

  return (
    <button
      type="button"
      className="rg-btn rg-btn--sm rg-digest__copy"
      onClick={copy}
      aria-label={label}
      title={label}
    >
      {copied ? 'Copied' : 'Copy'}
    </button>
  );
}
