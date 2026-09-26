/**
 * Resolves the bootstrap config once at startup and puts `{ config, client }`
 * in context. Children never see `undefined`: the provider renders a loading
 * state until the config is resolved.
 */

import { useEffect, useMemo, useState } from 'react';
import type { ReactNode } from 'react';
import {
  AppContext,
  createClient,
  resolveApiBase,
  resolveLocalApiBase,
} from '../services/clients';
import { loadUiConfig } from '../services/config';
import { LoadingState } from '../components/common/states';
import type { UiConfig } from '../types/ui';

export interface ConfigProviderProps {
  children: ReactNode;
  /** Tests (and the manage plan's stories) seed an explicit config. */
  config?: UiConfig;
}

export function ConfigProvider({ children, config: seeded }: ConfigProviderProps) {
  const [config, setConfig] = useState<UiConfig | undefined>(seeded);

  useEffect(() => {
    if (seeded) {
      setConfig(seeded);
      return;
    }
    let cancelled = false;
    // One attempt, then the read-only fallback. Never a retry loop.
    loadUiConfig().then((resolved) => {
      if (!cancelled) setConfig(resolved);
    });
    return () => {
      cancelled = true;
    };
  }, [seeded]);

  const value = useMemo(() => {
    if (!config) return undefined;
    return {
      config,
      client: createClient(resolveApiBase(config)),
      localClient: createClient(resolveLocalApiBase(config)),
    };
  }, [config]);

  if (!value) {
    return (
      <div className="p-8">
        <LoadingState label="Connecting to refgenie" />
      </div>
    );
  }

  return <AppContext.Provider value={value}>{children}</AppContext.Provider>;
}
