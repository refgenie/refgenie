import { Fragment } from 'react';
import { Link } from 'react-router-dom';

export interface Crumb {
  label: string;
  to?: string;
}

export interface BreadcrumbsProps {
  items: Crumb[];
}

export function Breadcrumbs({ items }: BreadcrumbsProps) {
  return (
    <nav aria-label="Breadcrumb" className="mb-4 text-sm">
      {items.map((item, index) => (
        <Fragment key={`${item.label}-${index}`}>
          {index > 0 && <span className="rg-muted"> / </span>}
          {item.to ? (
            <Link className="rg-link" to={item.to}>
              {item.label}
            </Link>
          ) : (
            <span className="rg-muted">{item.label}</span>
          )}
        </Fragment>
      ))}
    </nav>
  );
}
