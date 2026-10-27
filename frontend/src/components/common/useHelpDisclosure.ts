/**
 * The "?" toggle and its collapsible explainer, as props rather than markup.
 *
 * A component returning both pieces would not fit `FormField`, whose button
 * and hint live in different cells of a grid. A hook lets each caller place
 * them where its own layout needs while the accessibility wiring — the
 * `aria-controls`/`aria-expanded` pair, the button's name, and the fact that
 * collapsing is `sr-only` rather than removal — is written down once.
 *
 * Collapsed is a VISUAL state only. The hint stays in the DOM so it stays
 * readable to a screen reader and stays a valid `aria-describedby` target.
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
 * @param label Names the thing being explained, for the button's own name. A
 *   bare "Help" would make every toggle on the page indistinguishable in a
 *   screen reader's list of buttons.
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
