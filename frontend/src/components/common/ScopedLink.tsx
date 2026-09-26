/**
 * A `Link` that stays inside the current browse scope.
 *
 * Under `/local/*` the same browse components render against a connected local
 * refgenie. Their absolute links (`/genomes/{digest}`) would otherwise jump
 * back to the page's own backend carrying a digest that may not exist there.
 * Every genome/asset link inside a scopable component goes through this.
 *
 * Asset-class and recipe links stay plain `Link`s: there are no
 * `/local/asset-classes` or `/local/recipes` routes, so those links
 * intentionally exit the scope. If that ever proves confusing, the fix is to
 * add the routes, not to scope the links.
 */

import { Link } from 'react-router-dom';
import type { LinkProps } from 'react-router-dom';
import { useRouteBase } from '../../hooks/useRouteBase';

export function ScopedLink({ to, ...rest }: LinkProps) {
  const base = useRouteBase();
  const scoped = typeof to === 'string' && to.startsWith('/') ? `${base}${to}` : to;
  return <Link to={scoped} {...rest} />;
}
