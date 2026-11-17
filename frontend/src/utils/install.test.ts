// @vitest-environment node
import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { FALLBACK_VERSION, installCommand } from './install';

describe('installCommand', () => {
  it('pins the floor to the server version', () => {
    expect(installCommand('1.0.0a2')).toBe('pip install "refgenie>=1.0.0a2"');
  });

  it('adds an extra', () => {
    expect(installCommand('1.0.0a2', 'dash')).toBe('pip install "refgenie[dash]>=1.0.0a2"');
  });

  it('falls back when the version is unknown', () => {
    for (const version of ['unknown', '', null, undefined]) {
      expect(installCommand(version)).toBe(`pip install "refgenie>=${FALLBACK_VERSION}"`);
    }
  });

  it('drops a local build label', () => {
    expect(installCommand('1.0.0a3.dev4+g1a2b3c')).toBe('pip install "refgenie>=1.0.0a3.dev4"');
  });

  it('matches the install line in the no-JavaScript text of index.html', () => {
    const html = readFileSync(new URL('../../index.html', import.meta.url), 'utf8');
    expect(html).toContain(installCommand(FALLBACK_VERSION));
  });
});
