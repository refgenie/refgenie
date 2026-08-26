import { describe, expect, it, vi } from 'vitest';
import { fallbackConfig, loadUiConfig, parseServiceInfo } from './config';

const PAYLOAD = {
  id: 'org.refgenie.dash',
  name: 'refgenie',
  refgenie: {
    mode: 'local',
    api_base: '/v4',
    root_path: '',
    refgenie_version: '1.0.0a1',
    service_name: 'refgenie local dashboard',
    capabilities: { pull: true, remote_browse: true, downloads: false },
    links: { docs: 'https://docs.refgenie.org' },
  },
};

function response(body: unknown, status = 200): Response {
  return new Response(JSON.stringify(body), {
    status,
    headers: { 'content-type': 'application/json' },
  });
}

describe('parseServiceInfo', () => {
  it('reads the nested refgenie block', () => {
    const config = parseServiceInfo(PAYLOAD);
    expect(config.mode).toBe('local');
    expect(config.service_name).toBe('refgenie local dashboard');
    expect(config.degraded).toBe(false);
  });

  it('treats a missing capability key as false', () => {
    const config = parseServiceInfo(PAYLOAD);
    expect(config.capabilities.pull).toBe(true);
    expect(config.capabilities.jobs).toBe(false);
    expect(config.capabilities.archives).toBe(false);
  });

  it('merges default links so docs/github/openapi always resolve', () => {
    const config = parseServiceInfo(PAYLOAD);
    expect(config.links.openapi).toBe('/openapi.json');
  });
});

describe('loadUiConfig resolution chain', () => {
  it('uses the response when /service-info is found', async () => {
    const fetchImpl = vi.fn(async () => response(PAYLOAD));
    const config = await loadUiConfig(fetchImpl as never);
    expect(config.degraded).toBe(false);
    expect(config.mode).toBe('local');
  });

  it('falls back to read-only server mode on a 404', async () => {
    const fetchImpl = vi.fn(async () => response({ detail: 'Not Found' }, 404));
    const config = await loadUiConfig(fetchImpl as never);
    expect(config.degraded).toBe(true);
    expect(config.mode).toBe('server');
    expect(config.capabilities.pull).toBe(false);
    expect(config.capabilities.downloads).toBe(true);
    expect(fetchImpl).toHaveBeenCalledTimes(1);
  });

  it('falls back without crashing or retrying on malformed JSON', async () => {
    const fetchImpl = vi.fn(
      async () => new Response('not json', { headers: { 'content-type': 'application/json' } }),
    );
    const config = await loadUiConfig(fetchImpl as never);
    expect(config.degraded).toBe(true);
    expect(fetchImpl).toHaveBeenCalledTimes(1);
  });

  it('falls back when the payload has no refgenie block', async () => {
    const fetchImpl = vi.fn(async () => response({ id: 'x', name: 'y' }));
    const config = await loadUiConfig(fetchImpl as never);
    expect(config.degraded).toBe(true);
  });

  it('falls back on a network failure', async () => {
    const fetchImpl = vi.fn(async () => {
      throw new TypeError('Failed to fetch');
    });
    const config = await loadUiConfig(fetchImpl as never);
    expect(config.degraded).toBe(true);
  });
});

describe('VITE_API_BASE override', () => {
  it('wins over the api_base from /service-info', async () => {
    vi.stubEnv('VITE_API_BASE', 'https://api.refgenie.org/v4');
    try {
      const fetchImpl = vi.fn(async () => response(PAYLOAD));
      const config = await loadUiConfig(fetchImpl as never);
      expect(config.api_base).toBe('https://api.refgenie.org/v4');
    } finally {
      vi.unstubAllEnvs();
    }
  });

  it('is the fallback api_base when /service-info is unreachable', async () => {
    vi.stubEnv('VITE_API_BASE', 'https://api.refgenie.org/v4');
    try {
      expect(fallbackConfig().api_base).toBe('https://api.refgenie.org/v4');
    } finally {
      vi.unstubAllEnvs();
    }
  });
});

describe('fallbackConfig', () => {
  it('is read-only: every write capability is off', () => {
    const config = fallbackConfig();
    for (const key of [
      'pull',
      'build',
      'delete',
      'aliases_write',
      'subscriptions',
      'remote_browse',
      'genome_init',
      'jobs',
    ] as const) {
      expect(config.capabilities[key]).toBe(false);
    }
  });
});
