import { Icon } from './Icon';

export interface DownloadLinkProps {
  href: string;
  children: React.ReactNode;
  className?: string;
}

/**
 * An anchor whose click starts a byte transfer rather than a navigation.
 *
 * No `target="_blank"`: the server sends `Content-Disposition: attachment`, so
 * the current tab is never navigated away and a new tab would flash open and
 * close.
 *
 * No `download` attribute either. It is ignored cross-origin, which is exactly
 * the localhost-bridge case, so relying on it would be a lie on the one path
 * where it matters.
 *
 * The `sr-only` prefix, not a suffix: it makes the accessible name read
 * "Download hg38.fa.fai", which is the order a screen-reader user needs the
 * verb in.
 */
export function DownloadLink({ href, children, className = 'rg-link' }: DownloadLinkProps) {
  return (
    <a className={className} href={href}>
      <span className="sr-only">Download </span>
      {children}
      <Icon name="download" className="rg-icon--trailing" />
    </a>
  );
}
