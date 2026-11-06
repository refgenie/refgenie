import { useContext } from 'react';
import { AppContext } from '../services/clients';
import type { ApiClient } from '../services/http';

function useAppContext() {
  const value = useContext(AppContext);
  if (!value) {
    throw new Error('useApiClient must be used inside <ConfigProvider>');
  }
  return value;
}

/** The shared read API client (`/v4`). */
export function useApiClient(): ApiClient {
  return useAppContext().client;
}

/** The local-only surface client (`/v1`): remote browse, actions, jobs. */
export function useLocalApiClient(): ApiClient {
  return useAppContext().localClient;
}
