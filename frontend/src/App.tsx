import { useMemo } from 'react';
import { RouterProvider, createBrowserRouter } from 'react-router-dom';
import { createRoutes } from './app/routes';
import { useUiConfig } from './hooks/useUiConfig';

export function App() {
  const config = useUiConfig();
  // `root_path` becomes the router basename, for sub-path deployments.
  const router = useMemo(
    () =>
      createBrowserRouter(createRoutes(config), {
        basename: config.root_path || '/',
      }),
    [config],
  );
  return <RouterProvider router={router} />;
}
