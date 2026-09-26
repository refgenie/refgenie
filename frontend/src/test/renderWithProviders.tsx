import { render } from '@testing-library/react';
import { MemoryRouter } from 'react-router-dom';
import type { ReactElement, ReactNode } from 'react';
import { ConfigProvider } from '../app/ConfigProvider';
import { ResourceCacheContext } from '../services/resourceCacheContext';
import { createResourceCache, NO_RETRY } from '../services/resourceCache';
import { ToastProvider } from '../components/common/ToastProvider';
import { CAPABILITY_KEYS } from '../types/ui';
import type { Capabilities, UiConfig } from '../types/ui';

function capabilities(overrides: Partial<Capabilities> = {}): Capabilities {
  const base = Object.fromEntries(
    CAPABILITY_KEYS.map((key) => [key, false]),
  ) as unknown as Capabilities;
  return { ...base, ...overrides };
}

export const localConfig: UiConfig = {
  mode: 'local',
  api_base: '/v4',
  root_path: '',
  refgenie_version: '1.0.0a1',
  service_name: 'refgenie local dashboard',
  capabilities: capabilities({
    pull: true,
    build: true,
    delete: true,
    aliases_write: true,
    subscriptions: true,
    data_channels: true,
    remote_browse: true,
    genome_init: true,
    jobs: true,
    jobs_cancel: true,
  }),
  web_ui: null,
  degraded: false,
};

export const serverConfig: UiConfig = {
  ...localConfig,
  mode: 'server',
  service_name: 'refgenie server',
  capabilities: capabilities({ downloads: true, archives: true, seqcol: true, drs: true }),
};

export interface RenderOptions {
  config?: UiConfig;
  route?: string;
}

/**
 * The provider stack on its own, for `renderHook`.
 *
 * A fresh cache per render: the cache is a live object, so sharing one across
 * tests in a file would let a previous test's rows satisfy the next one's read.
 * Retries are off so a deliberate 500 fails immediately instead of backing off.
 */
export function createWrapper(options: RenderOptions = {}) {
  const cache = createResourceCache({ retry: NO_RETRY });

  return function Wrapper({ children }: { children: ReactNode }) {
    return (
      <ResourceCacheContext.Provider value={cache}>
        <ConfigProvider config={options.config ?? localConfig}>
          <ToastProvider>
            <MemoryRouter initialEntries={[options.route ?? '/']}>{children}</MemoryRouter>
          </ToastProvider>
        </ConfigProvider>
      </ResourceCacheContext.Provider>
    );
  };
}

export function renderWithProviders(ui: ReactElement, options: RenderOptions = {}) {
  return render(ui, { wrapper: createWrapper(options) });
}
