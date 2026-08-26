/**
 * Route table.
 *
 * Data loading stays in components via TanStack Query hooks — no router
 * `loader`s, because a loader cannot reach the React-context ApiClient.
 * The `*` catch-all is a CHILD of the layout so a 404 keeps nav and footer.
 *
 * The managed surfaces are registered only when their capability flags are on.
 * Nothing relies on CSS hiding: in server mode the actions and jobs routers are
 * never mounted, so a hand-typed URL has nothing behind it either.
 */

import type { RouteObject } from 'react-router-dom';
import { Navigate, useParams } from 'react-router-dom';
import { AppLayout } from '../components/layout/AppLayout';
import { RouteErrorPage } from '../pages/RouteErrorPage';
import { NotFoundPage } from '../pages/NotFoundPage';
import { GenomesPage } from '../pages/GenomesPage';
import { GenomePage } from '../pages/GenomePage';
import { AssetsPage } from '../pages/AssetsPage';
import { AssetPage } from '../pages/AssetPage';
import { AssetGroupPage } from '../pages/AssetGroupPage';
import { AssetClassesPage } from '../pages/AssetClassesPage';
import { AssetClassPage } from '../pages/AssetClassPage';
import { RecipesPage } from '../pages/RecipesPage';
import { RecipePage } from '../pages/RecipePage';
import { AliasesPage } from '../pages/AliasesPage';
import { RemotePage } from '../pages/RemotePage';
import { PullPage } from '../pages/PullPage';
import { AboutPage } from '../pages/AboutPage';
import { BuildPage } from '../pages/BuildPage';
import { ManagePage } from '../pages/ManagePage';
import { JobsPage } from '../pages/JobsPage';
import type { Capabilities } from '../types/ui';

/**
 * `/genomes/:digest/:assetDigest` is a documented bridge deep link (see the
 * contract next to `SPA_CLIENT_ROUTES` in `refgenie/server/const.py`). The
 * asset detail page is canonical at `/assets/:digest`, so forward there.
 */
function AssetDeepLink() {
  const { assetDigest } = useParams<{ assetDigest: string }>();
  return <Navigate to={`/assets/${assetDigest}`} replace />;
}

function managedRoutes(capabilities: Partial<Capabilities>): RouteObject[] {
  const routes: RouteObject[] = [];
  if (capabilities.build) {
    routes.push({ path: 'build', element: <BuildPage />, errorElement: <RouteErrorPage /> });
  }
  if (
    capabilities.subscriptions ||
    capabilities.aliases_write ||
    capabilities.genome_init ||
    capabilities.build
  ) {
    routes.push({ path: 'manage', element: <ManagePage />, errorElement: <RouteErrorPage /> });
  }
  if (capabilities.jobs) {
    routes.push({ path: 'jobs', element: <JobsPage />, errorElement: <RouteErrorPage /> });
  }
  return routes;
}

export function createRoutes(capabilities: Partial<Capabilities>): RouteObject[] {
  return [
    {
      path: '/',
      element: <AppLayout />,
      errorElement: <RouteErrorPage />,
      children: [
        // `/` renders the genome index; `/genomes` is the canonical list route.
        { index: true, element: <GenomesPage />, errorElement: <RouteErrorPage /> },
        { path: 'genomes', element: <GenomesPage />, errorElement: <RouteErrorPage /> },
        { path: 'genomes/:digest', element: <GenomePage />, errorElement: <RouteErrorPage /> },
        {
          path: 'genomes/:digest/:assetDigest',
          element: <AssetDeepLink />,
          errorElement: <RouteErrorPage />,
        },
        { path: 'assets', element: <AssetsPage />, errorElement: <RouteErrorPage /> },
        { path: 'assets/:digest', element: <AssetPage />, errorElement: <RouteErrorPage /> },
        { path: 'asset-groups/:id', element: <AssetGroupPage />, errorElement: <RouteErrorPage /> },
        { path: 'asset-classes', element: <AssetClassesPage />, errorElement: <RouteErrorPage /> },
        { path: 'asset-classes/:id', element: <AssetClassPage />, errorElement: <RouteErrorPage /> },
        { path: 'recipes', element: <RecipesPage />, errorElement: <RouteErrorPage /> },
        { path: 'recipes/:id', element: <RecipePage />, errorElement: <RouteErrorPage /> },
        { path: 'aliases', element: <AliasesPage />, errorElement: <RouteErrorPage /> },
        { path: 'remote', element: <RemotePage />, errorElement: <RouteErrorPage /> },
        // Registered unconditionally, like `remote`: a bridge hand-off URL
        // must degrade cleanly (in-page not-available), never 404.
        { path: 'pull', element: <PullPage />, errorElement: <RouteErrorPage /> },
        { path: 'about', element: <AboutPage />, errorElement: <RouteErrorPage /> },
        ...managedRoutes(capabilities),
        { path: '*', element: <NotFoundPage /> },
      ],
    },
  ];
}

/**
 * SPA client routes, kept in lockstep with `SPA_CLIENT_ROUTES` in
 * `refgenie/server/const.py`. The backend enforces both halves of the mirror:
 * `test_spa_route_mirror` (tests/test_app_modes.py) parses this const and the
 * registered route paths above and fails on any drift, and the non-collision
 * guard asserts the Python tuple's entries never shadow an API prefix.
 */
export const SPA_CLIENT_ROUTES = [
  '/genomes',
  '/assets',
  '/asset-groups',
  '/asset-classes',
  '/recipes',
  '/aliases',
  '/remote',
  '/manage',
  '/jobs',
  '/build',
  '/pull',
  '/about',
] as const;
