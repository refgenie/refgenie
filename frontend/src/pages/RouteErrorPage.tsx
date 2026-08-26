import { isRouteErrorResponse, useRouteError } from 'react-router-dom';

/**
 * Top-level boundary for render errors. Backend 500s arrive as JSON and are
 * handled by ErrorState; this catches everything React throws.
 */
export function RouteErrorPage() {
  const error = useRouteError();
  const message = isRouteErrorResponse(error)
    ? `${error.status} ${error.statusText}`
    : error instanceof Error
      ? error.message
      : 'Unknown error';

  return (
    <div className="rg-state rg-state--error" role="alert">
      <h1 className="text-2xl font-semibold mb-2">Something went wrong</h1>
      <p className="text-sm">{message}</p>
    </div>
  );
}
