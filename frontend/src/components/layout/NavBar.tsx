import { Link, NavLink } from 'react-router-dom';
import { NAV } from '../../app/navigation';
import { useUiConfig } from '../../hooks/useUiConfig';
import { ModeBadge } from './ModeBadge';
import { cn } from '../../utils/cn';

export function NavBar() {
  const config = useUiConfig();
  const entries = NAV.filter(
    (entry) => !entry.capability || config.capabilities[entry.capability] === true,
  );

  return (
    <nav className="rg-nav" aria-label="Main">
      <Link className="rg-nav__brand" to="/genomes">
        <img className="rg-nav__logo" src="refgenie_logo.svg" alt="" />
        <span>refgenie</span>
        <ModeBadge />
      </Link>
      <ul className="rg-nav__list">
        {entries.map((entry) => (
          <li key={entry.to}>
            <NavLink
              to={entry.to}
              className={({ isActive }) =>
                cn('rg-nav__item', isActive && 'rg-nav__item--active')
              }
            >
              {entry.label}
            </NavLink>
          </li>
        ))}
      </ul>
    </nav>
  );
}
