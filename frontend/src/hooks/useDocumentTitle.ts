/**
 * Sets the tab title per route, suffixed with `service_name` so a local dash and
 * the public server are easy to tell apart. Pass `undefined` while data is pending.
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
