/**
 * Props for a "?" toggle and its collapsible hint, so callers can place the two
 * pieces anywhere. Collapsing is visual only (`sr-only`): the hint stays readable
 * to screen readers and a valid `aria-describedby` target.
 */

import { useState } from 'react';
import { cn } from '../../utils/cn';

export interface HelpDisclosure {
  open: boolean;
  /** Spread onto a `<button>`. */
  toggleProps: {
    type: 'button';
    className: string;
    'aria-expanded': boolean;
    'aria-controls': string;
    'aria-label': string;
    onClick: () => void;
  };
  /** Spread onto the element holding the hint text. */
  hintProps: { id: string; className: string };
}

/**
 * @param hintId Also what the described control points `aria-describedby` at.
 * @param label Names the thing explained, so each toggle has a distinct accessible name.
 * @param hintClassName Extra classes for the hint element.
 */
export function useHelpDisclosure(
  hintId: string,
  label: string,
  hintClassName?: string,
): HelpDisclosure {
  const [open, setOpen] = useState(false);
  return {
    open,
    toggleProps: {
      type: 'button',
      className: 'rg-field__help',
      'aria-expanded': open,
      'aria-controls': hintId,
      'aria-label': `Help with ${label}`,
      onClick: () => setOpen((current) => !current),
    },
    hintProps: {
      id: hintId,
      className: cn('rg-field__hint', hintClassName, !open && 'sr-only'),
    },
  };
}
