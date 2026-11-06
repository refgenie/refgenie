import { describe, expect, it } from 'vitest';
import { vi } from 'vitest';

import { validatePing, SUPPORTED_BRIDGE_VERSIONS } from './contract';
import { isSelfLocal, isWebKitOnly } from './environment';
import { forget, recall, remember } from './persistence';
import { probeLocal } from './probe';

const validPing = {
  service: 'refgenie',
  bridge_version: 1,
  mode: 'local',
  refgenie_version: '1.0.0',
  api_version: 'v4',
  instance_id: '8f2c0000-0000-0000-0000-000000000000',
  instance_label: 'local refgenie',
  bridge_mode: 'read',
  action_header: 'X-Refgenie-Action',
  capabilities: { pull: true, archives: false },
  bridge: { actions_cross_origin: false },
  counts: { genomes: 12, assets: 47 },
};

describe('validatePing', () => {
  it('accepts a valid ping and normalizes booleans', () => {
    const result = validatePing(validPing);
    expect(result.ok).toBe(true);
    if (result.ok) {
      expect(result.ping.capabilities.pull).toBe(true);
      expect(result.ping.capabilities.archives).toBe(false);
      expect(result.ping.bridge.actions_cross_origin).toBe(false);
      expect(result.ping.counts).toEqual({ genomes: 12, assets: 47 });
    }
  });

  it('rejects a wrong service', () => {
    const result = validatePing({ ...validPing, service: 'jupyter' });
    expect(result).toEqual({ ok: false, reason: 'wrong-service' });
  });

  it('rejects an unsupported (newer) bridge version', () => {
    const unknownVersion = Math.max(...SUPPORTED_BRIDGE_VERSIONS) + 1;
    const result = validatePing({ ...validPing, bridge_version: unknownVersion });
    expect(result).toEqual({ ok: false, reason: 'unsupported-version' });
  });

  it.each([
    null,
    'a string',
    42,
    [],
    {},
    { service: 'refgenie' },
    { ...validPing, bridge_version: '1' },
    { ...validPing, instance_id: undefined },
    { ...validPing, capabilities: 'all' },
    { ...validPing, bridge: null },
  ])('rejects malformed payload %#', (payload) => {
    const result = validatePing(payload);
    expect(result.ok).toBe(false);
  });

  it('coerces non-boolean capability values to false, never truthy', () => {
    const result = validatePing({
      ...validPing,
      capabilities: { pull: 'yes', build: 1, delete: true },
    });
    expect(result.ok).toBe(true);
    if (result.ok) {
      expect(result.ping.capabilities.pull).toBe(false);
      expect(result.ping.capabilities.build).toBe(false);
      expect(result.ping.capabilities.delete).toBe(true);
    }
  });

  it('drops an undeclared capability key and keeps every declared one', () => {
    const result = validatePing({
      ...validPing,
      capabilities: { pull: true, teleport: true },
    });
    expect(result.ok).toBe(true);
    if (result.ok) {
      expect(result.ping.capabilities.pull).toBe(true);
      // Every key in the shared contract is present, and only those.
      expect(result.ping.capabilities).not.toHaveProperty('teleport');
      expect(result.ping.capabilities.jobs).toBe(false);
      expect(result.ping.capabilities.drs).toBe(false);
    }
  });

  it('accepts a payload with no `counts` key at all', () => {
    // `response_model_exclude_none=True` OMITS the field when the count
    // queries fail; it is never sent as null.
    const { counts: _counts, ...withoutCounts } = validPing;
    const result = validatePing(withoutCounts);
    expect(result.ok).toBe(true);
    if (result.ok) expect(result.ping.counts).toBeUndefined();
  });
});

describe('environment detectors', () => {
  it('flags Safari (WebKit vendor, no Chromium markers)', () => {
    expect(
      isWebKitOnly({
        vendor: 'Apple Computer, Inc.',
        userAgent:
          'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 Version/17.4 Safari/605.1.15',
      }),
    ).toBe(true);
  });

  it('does not flag Chrome (which also reports WebKit lineage)', () => {
    expect(
      isWebKitOnly({
        vendor: 'Google Inc.',
        userAgent: 'Mozilla/5.0 AppleWebKit/537.36 Chrome/125.0.0.0 Safari/537.36',
      }),
    ).toBe(false);
  });

  it('does not flag Chrome on iOS (Apple vendor but CriOS marker)', () => {
    expect(
      isWebKitOnly({
        vendor: 'Apple Computer, Inc.',
        userAgent: 'Mozilla/5.0 (iPhone) AppleWebKit/605.1.15 CriOS/125.0 Safari/604.1',
      }),
    ).toBe(false);
  });

  it('detects self-local for both loopback spellings, port-sensitive', () => {
    expect(isSelfLocal(8080, { origin: 'http://localhost:8080' })).toBe(true);
    expect(isSelfLocal(8080, { origin: 'http://127.0.0.1:8080' })).toBe(true);
    expect(isSelfLocal(8080, { origin: 'http://localhost:9000' })).toBe(false);
    expect(isSelfLocal(8080, { origin: 'https://refgenie.org' })).toBe(false);
  });
});

const fakeStorage = () => {
  const map = new Map<string, string>();
  return {
    getItem: (key: string) => map.get(key) ?? null,
    setItem: (key: string, value: string) => void map.set(key, value),
    removeItem: (key: string) => void map.delete(key),
  };
};

describe('persistence', () => {
  it('round-trips a remembered connection', () => {
    const storage = fakeStorage();
    const connection = {
      port: 8080,
      instanceId: 'abc',
      connectedAt: '2026-08-12T00:00:00Z',
    };
    remember(connection, storage);
    expect(recall(storage)).toEqual(connection);
  });

  it('forget clears the remembered connection', () => {
    const storage = fakeStorage();
    remember({ port: 8080, instanceId: 'abc', connectedAt: 'now' }, storage);
    forget(storage);
    expect(recall(storage)).toBeNull();
  });

  it('recall tolerates garbage and wrong shapes', () => {
    const storage = fakeStorage();
    storage.setItem('refgenie.localBridge', 'not json {');
    expect(recall(storage)).toBeNull();
    storage.setItem('refgenie.localBridge', JSON.stringify({ port: 'eight' }));
    expect(recall(storage)).toBeNull();
  });
});

describe('probeLocal outcome classification', () => {
  const jsonResponse = (body: unknown, status = 200) =>
    ({
      ok: status >= 200 && status < 300,
      status,
      json: async () => body,
    }) as Response;

  it('classifies a valid answer as connected', async () => {
    const fetchImpl = vi.fn(
      async (_input: RequestInfo | URL, _init?: RequestInit) => jsonResponse(validPing),
    );
    const result = await probeLocal(8080, fetchImpl as unknown as typeof fetch);
    expect(result.kind).toBe('connected');
    if (result.kind === 'connected') {
      expect(result.baseUrl).toBe('http://localhost:8080');
    }
    // The probe must stay a CORS *simple* request: no custom headers,
    // credentials omitted.
    const init = fetchImpl.mock.calls[0][1];
    expect(init?.credentials).toBe('omit');
    expect(init?.headers).toBeUndefined();
  });

  it('probes exactly one port and exactly one endpoint — never a scan', async () => {
    const fetchImpl = vi.fn(
      async (_input: RequestInfo | URL, _init?: RequestInit) => jsonResponse({}, 500),
    );
    await probeLocal(8080, fetchImpl as unknown as typeof fetch);
    expect(fetchImpl).toHaveBeenCalledTimes(1);
    expect(fetchImpl.mock.calls[0][0]).toBe('http://localhost:8080/ping');
  });

  it('classifies a 404 as unreachable', async () => {
    const fetchImpl = vi.fn(async () => jsonResponse({}, 404));
    const result = await probeLocal(8080, fetchImpl as unknown as typeof fetch);
    expect(result).toEqual({ kind: 'unreachable' });
  });

  it('classifies a 500 as unreachable', async () => {
    const fetchImpl = vi.fn(async () => jsonResponse({ oops: true }, 500));
    const result = await probeLocal(8080, fetchImpl as unknown as typeof fetch);
    expect(result).toEqual({ kind: 'unreachable' });
  });

  it('classifies a non-refgenie answer as unsupported/wrong-service', async () => {
    const fetchImpl = vi.fn(async () => jsonResponse({ hello: 'world' }));
    const result = await probeLocal(8080, fetchImpl as unknown as typeof fetch);
    expect(result).toEqual({ kind: 'unsupported', reason: 'wrong-service' });
  });

  it('classifies a newer bridge version as unsupported/unsupported-version', async () => {
    const fetchImpl = vi.fn(async () =>
      jsonResponse({ ...validPing, bridge_version: 99 }),
    );
    const result = await probeLocal(8080, fetchImpl as unknown as typeof fetch);
    expect(result).toEqual({ kind: 'unsupported', reason: 'unsupported-version' });
  });

  it('classifies a rejected fetch as unreachable (opaque by design)', async () => {
    const fetchImpl = vi.fn(async () => {
      throw new TypeError('Failed to fetch');
    });
    const result = await probeLocal(8080, fetchImpl as unknown as typeof fetch);
    expect(result).toEqual({ kind: 'unreachable' });
  });
});
