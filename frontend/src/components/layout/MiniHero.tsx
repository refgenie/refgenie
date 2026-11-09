/**
 * The standard page head: breadcrumb, title, optional explainer, actions.
 *
 * It also sets the document title, so a page must render MiniHero in all of
 * its states (loading and error too) or the tab keeps the previous title.
 * Pass `documentTitle={false}` while the record is pending.
 */

import type { ReactNode } from 'react';
import { useDocumentTitle } from '../../hooks/useDocumentTitle';

export interface MiniHeroProps {
  /** The <h1>, and the document title unless `documentTitle` overrides it. */
  title: string;
  /**
   * Document title when it differs from the heading. `false` means the record
   * is still loading, so the tab shows just the service name.
   */
  documentTitle?: string | false;
  /**
   * One or two plain sentences defining the page's subject. List pages define
   * their concept; detail pages need one only if no list page does.
   */
  lede?: ReactNode;
  /** A `<Breadcrumbs>` element, rendered above the title. */
  breadcrumbs?: ReactNode;
  /** Buttons, links and jump navs, right-aligned on the title's line. */
  actions?: ReactNode;
  /** Extra head content under the title line, e.g. a detail page's key facts. */
  children?: ReactNode;
}

export function MiniHero({
  title,
  documentTitle,
  lede,
  breadcrumbs,
  actions,
  children,
}: MiniHeroProps) {
  useDocumentTitle(documentTitle === false ? undefined : (documentTitle ?? title));

  return (
    <header className="rg-minihero">
      {breadcrumbs}
      <div className="rg-minihero__bar">
        <h1 className="rg-minihero__title">{title}</h1>
        {actions && <div className="rg-minihero__actions">{actions}</div>}
      </div>
      {lede && <p className="rg-minihero__lede">{lede}</p>}
      {children}
    </header>
  );
}
