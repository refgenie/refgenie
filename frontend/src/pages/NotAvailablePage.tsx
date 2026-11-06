import { useDocumentTitle } from '../hooks/useDocumentTitle';

export interface NotAvailablePageProps {
  title: string;
  reason: string;
}

/**
 * Shown for routes not enabled on this instance. Sets the document title
 * itself, since gated pages return this instead of their MiniHero.
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
