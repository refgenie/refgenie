/**
 * Route table. Data loads in components, not router loaders, because loaders
 * cannot reach the context ApiClient. Managed routes exist only when their
 * capability is on.
 */

import { lazy, Suspense, useContext } from 'react';
import type { RouteObject } from 'react-router-dom';
import { Navigate, Outlet, useParams } from 'react-router-dom';
import { AppLayout } from '../components/layout/AppLayout';
import { RouteErrorPage } from '../pages/RouteErrorPage';
import { NotFoundPage } from '../pages/NotFoundPage';
import { LandingPage } from '../pages/LandingPage';
import { GenomesPage } from '../pages/GenomesPage';
import { GenomePage } from '../pages/GenomePage';
import { SpeciesPage } from '../pages/SpeciesPage';
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
import { ServersPage } from '../pages/ServersPage';
import { JobsPage } from '../pages/JobsPage';
import { CollectionsPage } from '../pages/CollectionsPage';
import { CollectionPage } from '../pages/CollectionPage';
import { LoadingState } from '../components/common/states';
import { LocalScope } from '../components/bridge/LocalScope';
import { BridgeScopeContext } from '../services/bridge/scope';
import type { Capabilities, UiConfig } from '../types/ui';

// Lazy: d3 and the bundled taxonomy (~200 KB) are only needed on `/tree`.
const TreeExplorerPage = lazy(() =>
  import('../pages/TreeExplorerPage').then((m) => ({ default: m.TreeExplorerPage })),
);

/** Redirects the bridge deep link `/genomes/:digest/:assetDigest` to `/assets/:digest`. */
function AssetDeepLink() {
  const { assetDigest } = useParams<{ assetDigest: string }>();
  const { routeBase } = useContext(BridgeScopeContext);
  return <Navigate to={`${routeBase}/assets/${assetDigest}`} replace />;
}

/**
 * One genome page per digest. Keyed, so moving to another genome remounts it:
 * no download, extract or compare state (or job in flight) carries over.
 */
export function GenomeRoute() {
  const { digest } = useParams<{ digest: string }>();
  return <GenomePage key={digest} />;
}

/** Browse pages, mounted at the root and again under `/local`. */
const BROWSE_ROUTES: RouteObject[] = [
  { path: 'genomes', element: <GenomesPage /> },
  { path: 'genomes/:digest', element: <GenomeRoute /> },
  { path: 'genomes/:digest/:assetDigest', element: <AssetDeepLink /> },
  { path: 'species', element: <SpeciesPage /> },
  { path: 'assets/:digest', element: <AssetPage /> },
  { path: 'asset-groups/:id', element: <AssetGroupPage /> },
  { path: 'asset-classes', element: <AssetClassesPage /> },
  { path: 'asset-classes/:id', element: <AssetClassPage /> },
  { path: 'recipes', element: <RecipesPage /> },
  { path: 'recipes/:id', element: <RecipePage /> },
  { path: 'aliases', element: <AliasesPage /> },
];

/**
 * `/local/*`: read-only browse pages against a connected local refgenie.
 * Skipped on a local dash, which already is that refgenie. `LocalScope` masks
 * write capabilities.
 */
function localRoutes(config: UiConfig): RouteObject[] {
  if (config.mode === 'local') return [];
  return [
    {
      path: 'local',
      element: (
        <LocalScope>
          <Outlet />
        </LocalScope>
      ),
      children: [
        { index: true, element: <GenomesPage /> },
        ...BROWSE_ROUTES,
        // Read-only here; subscribing happens in the local dash.
        { path: 'servers', element: <ServersPage /> },
        { path: 'jobs', element: <JobsPage /> },
      ],
    },
  ];
}

function managedRoutes(capabilities: Partial<Capabilities>): RouteObject[] {
  const routes: RouteObject[] = [];
  if (capabilities.build) routes.push({ path: 'build', element: <BuildPage /> });
  if (capabilities.subscriptions || capabilities.remote_browse) {
    routes.push({ path: 'servers', element: <ServersPage /> });
  }
  if (
    capabilities.subscriptions ||
    capabilities.aliases_write ||
    capabilities.genome_init ||
    capabilities.build
  ) {
    routes.push({ path: 'manage', element: <ManagePage /> });
  }
  if (capabilities.jobs) routes.push({ path: 'jobs', element: <JobsPage /> });
  return routes;
}

export function createRoutes(config: UiConfig): RouteObject[] {
  const capabilities: Partial<Capabilities> = config.capabilities;
  return [
    {
      path: '/',
      element: <AppLayout />,
      errorElement: <RouteErrorPage />,
      children: [
        {
          // Pathless, so a page error renders inside the layout.
          errorElement: <RouteErrorPage />,
          children: [
            { index: true, element: <LandingPage /> },
            ...BROWSE_ROUTES,
            {
              // Always registered so a shared link resolves; only the nav entry is gated.
              path: 'tree',
              element: (
                <Suspense fallback={<LoadingState label="Loading the tree" />}>
                  <TreeExplorerPage />
                </Suspense>
              ),
            },
            // Server collections; not under /local. Always registered so a shared link resolves.
            { path: 'collections', element: <CollectionsPage /> },
            { path: 'collections/:id', element: <CollectionPage /> },
            { path: 'remote', element: <RemotePage /> },
            // Always registered: a bridge hand-off URL shows "not available", never a 404.
            { path: 'pull', element: <PullPage /> },
            { path: 'about', element: <AboutPage /> },
            ...managedRoutes(capabilities),
            ...localRoutes(config),
            { path: '*', element: <NotFoundPage /> },
          ],
        },
      ],
    },
  ];
}

/**
 * Must match `SPA_CLIENT_ROUTES` in `refgenie/server/const.py`;
 * `test_spa_route_mirror` checks this.
 */
export const SPA_CLIENT_ROUTES = [
  '/genomes',
  '/species',
  '/assets',
  '/asset-groups',
  '/asset-classes',
  '/recipes',
  '/aliases',
  '/tree',
  '/collections',
  '/local',
  '/remote',
  '/manage',
  '/servers',
  '/jobs',
  '/build',
  '/pull',
  '/about',
] as const;
