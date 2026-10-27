/**
 * The standard page head: breadcrumb, title, optional one-or-two-sentence
 * explainer, optional actions.
 *
 * The landing page's hero at inner-page scale -- see the block comment on
 * `.rg-minihero` in components.css for why it is a header rather than a banner.
 *
 * It owns the document title, because the document title and the <h1> are the
 * same fact and should not be two independent statements in every page file.
 * That has one consequence worth knowing: a page that early-returns a loading
 * or error state before rendering MiniHero sets no document title at all. Every
 * detail page therefore renders MiniHero in ALL of its states and passes
 * `documentTitle={false}` while its record is still pending -- which shows the
 * bare service name, never the previous page's subject.
 *
 * The explainer rule the `lede` prop exists for: a LIST page defines the
 * concept it lists, once. A DETAIL page carries an explainer only when no list
 * page above it defines the word -- which today is `/assets/:digest` and
 * `/asset-groups/:id` and nothing else.
 */

import type { ReactNode } from 'react';
import { useDocumentTitle } from '../../hooks/useDocumentTitle';

export interface MiniHeroProps {
  /** The <h1>, and the document title unless `documentTitle` overrides it. */
  title: string;
  /**
   * Document title, when it legitimately differs from the heading -- `Remote`
   * for a page headed "Remote assets", or a group name qualified by its genome.
   * `false` means the caller's record is still loading: fall back to the bare
   * service name rather than titling the tab with a placeholder.
   */
  documentTitle?: string | false;
  /**
   * One or two plain sentences saying what this page's subject IS, written for
   * someone who has never heard the words "asset class", "asset group",
   * "recipe" or "digest".
   */
  lede?: ReactNode;
  /** A `<Breadcrumbs>` element. Rendered above the title, gapped by the block. */
  breadcrumbs?: ReactNode;
  /** Buttons, links and jump navs, right-aligned on the title's line. */
  actions?: ReactNode;
}

export function MiniHero({ title, documentTitle, lede, breadcrumbs, actions }: MiniHeroProps) {
  useDocumentTitle(documentTitle === false ? undefined : (documentTitle ?? title));

  return (
    <header className="rg-minihero">
      {breadcrumbs}
      <div className="rg-minihero__bar">
        <h1 className="rg-minihero__title">{title}</h1>
        {actions && <div className="rg-minihero__actions">{actions}</div>}
      </div>
      {lede && <p className="rg-minihero__lede">{lede}</p>}
    </header>
  );
}
