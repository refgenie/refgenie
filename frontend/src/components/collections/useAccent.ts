import { useCallback } from 'react';

/**
 * A ref callback that sets `--rg-collection-accent` on its element. The
 * color is server config, so it cannot be a class; the lint rule bans the
 * `style` attribute, so it is set through the DOM instead.
 */
export function useAccent<T extends HTMLElement>(color: string) {
  return useCallback(
    (element: T | null) => {
      element?.style.setProperty('--rg-collection-accent', color);
    },
    [color],
  );
}
