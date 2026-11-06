import { Icon } from './Icon';

export interface DownloadLinkProps {
  href: string;
  children: React.ReactNode;
  className?: string;
}

/**
 * A link that downloads a file. It needs no `target` or `download` attribute:
 * the server sends `Content-Disposition: attachment`, and `download` is ignored
 * cross-origin anyway. The `sr-only` prefix reads as "Download <file>".
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
