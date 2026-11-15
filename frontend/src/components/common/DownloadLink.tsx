import { Icon } from './Icon';
import { saveUrlAs } from '../../services/saveStream';

export interface DownloadLinkProps {
  href: string;
  children: React.ReactNode;
  className?: string;
  /** Save the file under this name. The page reads it into memory: small files only. */
  saveAs?: string;
}

/**
 * A link that downloads a file. It needs no `target` or `download` attribute:
 * the server sends `Content-Disposition: attachment`, and `download` is ignored
 * cross-origin anyway. The `sr-only` prefix reads as "Download <file>".
 */
export function DownloadLink({
  href,
  children,
  className = 'rg-link',
  saveAs,
}: DownloadLinkProps) {
  const save = (event: React.MouseEvent) => {
    if (!saveAs) return;
    event.preventDefault();
    // If the page cannot fetch it, the plain link still downloads it.
    saveUrlAs(href, saveAs).catch(() => window.location.assign(href));
  };
  return (
    <a className={className} href={href} onClick={save}>
      <span className="sr-only">Download </span>
      {children}
      <Icon name="download" className="rg-icon--trailing" />
    </a>
  );
}
