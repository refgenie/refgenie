import { describe, expect, it } from 'vitest';
import type { RouteObject } from 'react-router-dom';
import { createRoutes } from './routes';
import { localConfig, serverConfig } from '../test/renderWithProviders';

function paths(routes: RouteObject[], prefix = ''): string[] {
  return routes.flatMap((route) => {
    const path =
      route.path && route.path !== '/'
        ? `${prefix}${route.path}`
        : prefix.replace(/\/$/, '');
    return [
      ...(path ? [path] : []),
      ...(route.children ? paths(route.children, path ? `${path}/` : '') : []),
    ];
  });
}

describe('createRoutes', () => {
  it('registers the collection pages at the root, not under /local', () => {
    const registered = paths(createRoutes(serverConfig));
    expect(registered).toContain('collections');
    expect(registered).toContain('collections/:id');
    expect(registered).not.toContain('local/collections');
  });

  it('registers the /local branch on a page that is not itself the dash', () => {
    const registered = paths(createRoutes(serverConfig));
    expect(registered).toContain('local');
    expect(registered).toContain('local/genomes/:digest');
    expect(registered).toContain('local/assets/:digest');
  });

  it('registers the read-only /local views for classes, recipes, aliases and jobs', () => {
    const registered = paths(createRoutes(serverConfig));
    for (const path of [
      'local/asset-classes',
      'local/asset-classes/:id',
      'local/recipes',
      'local/recipes/:id',
      'local/aliases',
      'local/jobs',
    ]) {
      expect(registered).toContain(path);
    }
  });

  it('registers no /local route on a local dash — it IS the local refgenie', () => {
    const registered = paths(createRoutes(localConfig));
    expect(registered.filter((path) => path.startsWith('local'))).toEqual([]);
  });

  it('registers /species in both modes', () => {
    expect(paths(createRoutes(serverConfig))).toContain('species');
    expect(paths(createRoutes(localConfig))).toContain('species');
  });

  it('registers the /local species twin only off the dash', () => {
    expect(paths(createRoutes(serverConfig))).toContain('local/species');
    expect(paths(createRoutes(localConfig))).not.toContain('local/species');
  });
});
