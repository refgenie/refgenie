import { Link } from 'react-router-dom';

export function NotFoundPage() {
  return (
    <div className="rg-state rg-state--empty">
      <h1 className="text-2xl font-semibold mb-2">Page not found</h1>
      <p className="mb-4">There is no page at this address.</p>
      <Link className="rg-link" to="/genomes">
        Back to genomes
      </Link>
    </div>
  );
}
