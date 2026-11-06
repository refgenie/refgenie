import { Fragment } from 'react';
import { ScopedLink } from './ScopedLink';

export interface Crumb {
  label: string;
  to?: string;
}

export interface BreadcrumbsProps {
  items: Crumb[];
}

export function Breadcrumbs({ items }: BreadcrumbsProps) {
  return (
    <nav aria-label="Breadcrumb" className="text-sm">
      {items.map((item, index) => (
        <Fragment key={`${item.label}-${index}`}>
          {index > 0 && <span className="rg-muted"> / </span>}
          {item.to ? (
            <ScopedLink className="rg-link" to={item.to}>
              {item.label}
            </ScopedLink>
          ) : (
            <span className="rg-muted">{item.label}</span>
          )}
        </Fragment>
      ))}
    </nav>
  );
}
