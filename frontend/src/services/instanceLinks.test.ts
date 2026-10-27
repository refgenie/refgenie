import { describe, expect, it } from 'vitest';
import { absoluteServerRoot, instanceUrl, serverRoot } from './instanceLinks';

describe('serverRoot', () => {
  it('drops the version suffix from a cross-origin base', () => {
    expect(serverRoot('https://api.refgenie.org/v4')).toBe('https://api.refgenie.org');
  });

  it('empties a same-origin base, so the URL stays root-relative', () => {
    expect(serverRoot('/v4')).toBe('');
  });

  it('keeps a sub-path root, which resolveApiBase already folded in', () => {
    expect(serverRoot('/refgenie/v4')).toBe('/refgenie');
  });
});

describe('absoluteServerRoot', () => {
  it('names the page origin when the API is same-origin', () => {
    // The bridge sends this to another process: `''` would tell the local
    // refgenie nothing, and it would fall back to its own subscriptions.
    expect(absoluteServerRoot('/v4', 'https://refgenie.org')).toBe('https://refgenie.org');
  });

  it('keeps a cross-origin API base as it is', () => {
    expect(absoluteServerRoot('https://api.refgenie.org/v4', 'https://refgenie.org')).toBe(
      'https://api.refgenie.org',
    );
  });

  it('carries a sub-path deployment onto the origin', () => {
    expect(absoluteServerRoot('/refgenie/v4', 'https://example.org')).toBe(
      'https://example.org/refgenie',
    );
  });
});

describe('instanceUrl', () => {
  it('reaches the API origin, not the SPA origin', () => {
    expect(instanceUrl('https://api.refgenie.org/v4', '/openapi.json')).toBe(
      'https://api.refgenie.org/openapi.json',
    );
  });

  it('stays relative when the API is same-origin', () => {
    expect(instanceUrl('/v4', '/openapi.json')).toBe('/openapi.json');
  });

  it('carries the sub-path', () => {
    expect(instanceUrl('/refgenie/v4', '/docs')).toBe('/refgenie/docs');
  });
});
