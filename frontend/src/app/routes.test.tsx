import { describe, expect, it } from 'vitest';
import type { RouteObject } from 'react-router-dom';
import { SPA_CLIENT_ROUTES, createRoutes } from './routes';
import { localConfig, serverConfig } from '../test/renderWithProviders';

function paths(routes: RouteObject[]): string[] {
  return routes.flatMap((route) => [
    ...(route.path ? [route.path] : []),
    ...(route.children ? paths(route.children) : []),
  ]);
}

describe('createRoutes', () => {
  it('registers the /local branch on a page that is not itself the dash', () => {
    const registered = paths(createRoutes(serverConfig));
    expect(registered).toContain('local');
    expect(registered).toContain('local/genomes/:digest');
    expect(registered).toContain('local/assets/:digest');
  });

  it('registers no /local route on a local dash — it IS the local refgenie', () => {
    const registered = paths(createRoutes(localConfig));
    expect(registered.filter((path) => path.startsWith('local'))).toEqual([]);
  });

  it('registers /species in both modes and mirrors it', () => {
    expect(paths(createRoutes(serverConfig))).toContain('species');
    expect(paths(createRoutes(localConfig))).toContain('species');
    expect(SPA_CLIENT_ROUTES).toContain('/species');
  });

  it('registers the /local species twin only off the dash', () => {
    expect(paths(createRoutes(serverConfig))).toContain('local/species');
    expect(paths(createRoutes(localConfig))).not.toContain('local/species');
  });

  it('keeps /local in the catch-all mirror', () => {
    // `tests/server/test_app.py::test_spa_route_mirror` asserts the other half of
    // this against `refgenie/server/const.py`.
    expect(SPA_CLIENT_ROUTES).toContain('/local');
  });
});
