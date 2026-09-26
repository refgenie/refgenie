import react from '@vitejs/plugin-react';
import { loadEnv } from 'vite';
import { defineConfig } from 'vitest/config';

/**
 * Paths that belong to the backend, never to the SPA router. In dev these are
 * proxied to a running `refgenie dash` / `refgenie serve`; in production they
 * are the prefixes the server's SPA catch-all refuses to swallow.
 */
const BACKEND_PATHS = [
  '/v4',
  '/v1',
  '/service-info',
  '/ping',
  '/openapi.json',
  '/seqcol',
  '/ga4gh',
  '/data_channel',
];

export default defineConfig(({ mode }) => {
  const env = loadEnv(mode, process.cwd(), '');
  // 8080 is the `refgenie dash` default (cli/commands/serve.py::DashModel).
  // Set VITE_DEV_BACKEND=http://127.0.0.1:8000 to develop against `refgenie serve`.
  const target = env.VITE_DEV_BACKEND ?? 'http://127.0.0.1:8080';

  return {
    plugins: [react()],
    // Absolute, not './', so deep routes like /genomes/<digest> resolve bundle
    // URLs correctly. Sub-path deployments are handled by spa.py rewriting the
    // literal <base href="/"> in index.html at app construction.
    base: '/',
    server: {
      port: 5173,
      proxy: Object.fromEntries(
        BACKEND_PATHS.map((p) => [p, { target, changeOrigin: true }]),
      ),
    },
    build: {
      outDir: '../refgenie/server/webui',
      emptyOutDir: true,
      // NOT the Vite default 'assets': that would put the hashed bundle at
      // /assets/*, colliding with the SPA's own /assets browse route and with
      // the /v4/assets API family.
      assetsDir: '_app',
      target: 'esnext',
      // Source maps must never ride into every `pip install`.
      sourcemap: false,
    },
    test: {
      globals: true,
      environment: 'jsdom',
      setupFiles: ['./src/test/setup.ts'],
      css: false,
      include: ['src/**/*.test.{ts,tsx}'],
    },
  };
});
