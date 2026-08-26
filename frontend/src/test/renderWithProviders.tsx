import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { render } from '@testing-library/react';
import { MemoryRouter } from 'react-router-dom';
import type { ReactElement, ReactNode } from 'react';
import { ConfigProvider } from '../app/ConfigProvider';
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
    remote_browse: true,
    genome_init: true,
    jobs: true,
    jobs_cancel: true,
  }),
  links: { docs: 'https://docs.refgenie.org', openapi: '/openapi.json' },
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

/** The provider stack on its own, for `renderHook`. */
export function createWrapper(options: RenderOptions = {}) {
  const client = new QueryClient({
    defaultOptions: { queries: { retry: false, gcTime: 0 } },
  });

  return function Wrapper({ children }: { children: ReactNode }) {
    return (
      <QueryClientProvider client={client}>
        <ConfigProvider config={options.config ?? localConfig}>
          <ToastProvider>
            <MemoryRouter initialEntries={[options.route ?? '/']}>{children}</MemoryRouter>
          </ToastProvider>
        </ConfigProvider>
      </QueryClientProvider>
    );
  };
}

export function renderWithProviders(ui: ReactElement, options: RenderOptions = {}) {
  return render(ui, { wrapper: createWrapper(options) });
}
