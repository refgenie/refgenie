import { useMemo } from 'react';
import { RouterProvider, createBrowserRouter } from 'react-router-dom';
import { createRoutes } from './app/routes';
import { useUiConfig } from './hooks/useUiConfig';

export function App() {
  const config = useUiConfig();
  // `root_path` from /service-info becomes the router basename, which is what
  // makes a sub-path deployment work. The managed routes exist only when their
  // capability flags are on, so a server-mode instance has nothing behind them.
  const router = useMemo(
    () =>
      createBrowserRouter(createRoutes(config.capabilities), {
        basename: config.root_path || '/',
      }),
    [config.root_path, config.capabilities],
  );
  return <RouterProvider router={router} />;
}
