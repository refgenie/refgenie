/**
 * The tree's SVG text and stroke values, read from the design tokens.
 *
 * SVG `font-size` is an attribute set by d3 inside a zoom transform, so it
 * cannot be a CSS class. Resolving the tokens here is what keeps the tree on
 * the same type scale as the rest of the app instead of inventing px values.
 * Fallbacks match tokens.css and cover jsdom, where no stylesheet is applied.
 */

export interface TreeTheme {
  fontFamily: string;
  speciesSize: number;
  groupSize: number;
  groupMutedSize: number;
  speciesWeight: string;
  groupWeight: string;
  groupMutedWeight: string;
  text: string;
  textMuted: string;
  accent: string;
}

/** Default root font size when `getComputedStyle` gives nothing usable. */
const FALLBACK_ROOT_PX = 16;

function px(
  styles: CSSStyleDeclaration,
  token: string,
  fallbackRem: number,
  rootPx: number,
): number {
  const raw = styles.getPropertyValue(token).trim();
  const rem = raw.endsWith('rem') ? Number.parseFloat(raw) : NaN;
  return (Number.isFinite(rem) ? rem : fallbackRem) * rootPx;
}

function rgb(styles: CSSStyleDeclaration, token: string, fallback: string): string {
  const triplet = styles.getPropertyValue(token).trim();
  return triplet ? `rgb(${triplet})` : fallback;
}

function text(styles: CSSStyleDeclaration, token: string, fallback: string): string {
  return styles.getPropertyValue(token).trim() || fallback;
}

export function readTreeTheme(element: Element): TreeTheme {
  const styles = getComputedStyle(element);
  const rootSize = Number.parseFloat(getComputedStyle(document.documentElement).fontSize);
  const rootPx = Number.isFinite(rootSize) ? rootSize : FALLBACK_ROOT_PX;

  return {
    fontFamily: text(styles, '--font-family-base', 'system-ui, sans-serif'),
    speciesSize: px(styles, '--font-size-sm', 0.875, rootPx),
    groupSize: px(styles, '--font-size-base', 1, rootPx),
    groupMutedSize: px(styles, '--font-size-sm', 0.875, rootPx),
    speciesWeight: text(styles, '--font-weight-semibold', '600'),
    groupWeight: text(styles, '--font-weight-bold', '700'),
    groupMutedWeight: text(styles, '--font-weight-normal', '400'),
    text: rgb(styles, '--color-text', 'rgb(17, 24, 39)'),
    textMuted: rgb(styles, '--color-border-strong', 'rgb(173, 181, 189)'),
    accent: rgb(styles, '--color-accent', 'rgb(13, 110, 253)'),
  };
}
