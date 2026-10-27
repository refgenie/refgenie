/**
 * Per-route document title.
 *
 * `index.html` ships a single `<title>refgenie</title>`; without this hook
 * every tab, history entry and bookmark in the app would be indistinguishable
 * from every other one.
 *
 * The suffix is `service_name` from `/service-info` rather than a literal:
 * a local dash and the public server are then telling apart in a tab strip,
 * which is the whole point of having a title at all.
 *
 * Pages pass `undefined` while their own data is pending so a title never
 * flashes the previous page's subject.
 */

import { useEffect } from 'react';
import { useUiConfig } from './useUiConfig';

export function useDocumentTitle(title: string | undefined) {
  const { service_name } = useUiConfig();

  useEffect(() => {
    const previous = document.title;
    document.title = title ? `${title} · ${service_name}` : service_name;
    return () => {
      document.title = previous;
    };
  }, [title, service_name]);
}
