export interface NotAvailablePageProps {
  title: string;
  reason: string;
}

/**
 * Rendered for routes that exist but are not enabled on this instance, and for
 * the `/manage` and `/jobs` namespaces reserved by the frontend-manage plan.
 */
export function NotAvailablePage({ title, reason }: NotAvailablePageProps) {
  return (
    <div className="rg-state rg-state--empty">
      <h1 className="text-2xl font-semibold mb-2">{title}</h1>
      <p>{reason}</p>
    </div>
  );
}
