import { useDocumentTitle } from '../hooks/useDocumentTitle';

export interface NotAvailablePageProps {
  title: string;
  reason: string;
}

/**
 * Rendered for routes that exist but are not enabled on this instance, and for
 * the `/manage` and `/jobs` namespaces reserved by the frontend-manage plan.
 *
 * It sets the document title itself. Every capability-gated page early-returns
 * this instead of its `MiniHero`, and `MiniHero` is what owns the title
 * everywhere else, so without this call a gated route would leave the tab
 * reading the bare service name. This is deliberately NOT a mini-hero: it is an
 * `rg-state--empty` block with a heading inside, and that shape is right for a
 * page with no subject.
 */
export function NotAvailablePage({ title, reason }: NotAvailablePageProps) {
  useDocumentTitle(title);

  return (
    <div className="rg-state rg-state--empty">
      <h1 className="text-2xl font-semibold mb-2">{title}</h1>
      <p>{reason}</p>
    </div>
  );
}
