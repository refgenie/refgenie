/**
 * The project's entire icon set.
 *
 * Inline SVG, no dependency, no icon font: the local dash must render with no
 * network at all, and `style-guard` bans both Bootstrap and Tailwind assets.
 * Geometry is Feather's (MIT), redrawn on a 24x24 grid, stroke-only so that one
 * `currentColor` inherits every colour the icon will ever need.
 *
 * `aria-hidden` is hard-coded and has no prop to override it. An icon is never
 * the accessible name of anything; the fact a glyph carries reaches the
 * accessibility tree as text, from the component that owns the glyph. See
 * `ExternalLink` and `DownloadLink` for the two ways that is done.
 *
 * The set is deliberately small: the link audit produced exactly three link
 * behaviours, plus three job kinds and a confirmation check. If it ever passes
 * roughly 25 icons, that is the point to reopen the npm-dependency question.
 */

import { cn } from '../../utils/cn';

export type IconName =
  | 'check'
  | 'code'
  | 'copy'
  | 'download'
  | 'external'
  | 'plus'
  | 'wrench';

const PATHS: Record<IconName, React.ReactNode> = {
  check: <polyline points="20 6 9 17 4 12" />,
  code: (
    <>
      <polyline points="16 18 22 12 16 6" />
      <polyline points="8 6 2 12 8 18" />
    </>
  ),
  copy: (
    <>
      <rect x="9" y="9" width="13" height="13" rx="2" ry="2" />
      <path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1" />
    </>
  ),
  download: (
    <>
      <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4" />
      <polyline points="7 10 12 15 17 10" />
      <line x1="12" y1="15" x2="12" y2="3" />
    </>
  ),
  external: (
    <>
      <path d="M18 13v6a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V8a2 2 0 0 1 2-2h6" />
      <polyline points="15 3 21 3 21 9" />
      <line x1="10" y1="14" x2="21" y2="3" />
    </>
  ),
  plus: (
    <>
      <line x1="12" y1="5" x2="12" y2="19" />
      <line x1="5" y1="12" x2="19" y2="12" />
    </>
  ),
  wrench: (
    <path d="M14.7 6.3a1 1 0 0 0 0 1.4l1.6 1.6a1 1 0 0 0 1.4 0l3.77-3.77a6 6 0 0 1-7.94 7.94l-6.91 6.91a2.12 2.12 0 0 1-3-3l6.91-6.91a6 6 0 0 1 7.94-7.94l-3.76 3.76z" />
  ),
};

/** Every name in `IconName`, for tests and for iteration. */
export const ICON_NAMES = Object.keys(PATHS) as IconName[];

export interface IconProps {
  name: IconName;
  /** `rg-icon--trailing`, `rg-icon--lg`. Never a colour or a size. */
  className?: string;
}

export function Icon({ name, className }: IconProps) {
  return (
    <svg
      className={cn('rg-icon', className)}
      viewBox="0 0 24 24"
      fill="none"
      stroke="currentColor"
      strokeWidth={2}
      strokeLinecap="round"
      strokeLinejoin="round"
      aria-hidden="true"
      focusable="false"
    >
      {PATHS[name]}
    </svg>
  );
}
