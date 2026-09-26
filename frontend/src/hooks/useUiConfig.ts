import { useContext } from 'react';
import { AppContext } from '../services/clients';
import type { UiConfig } from '../types/ui';

/** The bootstrap config. Never undefined: the provider blocks until resolved. */
export function useUiConfig(): UiConfig {
  const value = useContext(AppContext);
  if (!value) {
    throw new Error('useUiConfig must be used inside <ConfigProvider>');
  }
  return value.config;
}
