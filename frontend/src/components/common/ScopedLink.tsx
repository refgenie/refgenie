/**
 * A `Link` that stays inside the current browse scope, so links under
 * `/local/*` do not jump back to the page's own backend.
 */

import { Link } from 'react-router-dom';
import type { LinkProps } from 'react-router-dom';
import { useRouteBase } from '../../hooks/useRouteBase';

export function ScopedLink({ to, ...rest }: LinkProps) {
  const base = useRouteBase();
  const scoped = typeof to === 'string' && to.startsWith('/') ? `${base}${to}` : to;
  return <Link to={scoped} {...rest} />;
}
