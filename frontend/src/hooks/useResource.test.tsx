/**
 * `useResource` at the React level: one cache entry shared by several
 * components, gating, refetch and cancellation on unmount.
 */

import { describe, expect, it, vi } from 'vitest';
import { useState } from 'react';
import { render, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import type { ReactNode } from 'react';
import { useResource } from './useResource';
import { createResourceCache } from '../services/resourceCache';
import { ResourceCacheContext } from '../services/resourceCacheContext';
import type { ResourceCache, ResourceFetcher } from '../services/resourceCache';

function wrapperFor(cache: ResourceCache) {
  return function Wrapper({ children }: { children: ReactNode }) {
    return (
      <ResourceCacheContext.Provider value={cache}>{children}</ResourceCacheContext.Provider>
    );
  };
}

interface ProbeProps {
  label: string;
  fetcher: ResourceFetcher<string>;
  params?: object;
  enabled?: boolean;
}

function Probe({ label, fetcher, params = {}, enabled = true }: ProbeProps) {
  const query = useResource(['genomes', params], fetcher, { enabled, staleTime: 60_000 });
  return (
    <div>
      <span data-testid={label}>{query.isPending ? 'pending' : String(query.data)}</span>
      <button type="button" onClick={() => void query.refetch()}>
        refresh {label}
      </button>
    </div>
  );
}

describe('useResource', () => {
  it('serves two components from one request', async () => {
    const cache = createResourceCache();
    const fetcher = vi.fn(async () => 'rows');

    render(
      <>
        <Probe label="a" fetcher={fetcher} />
        <Probe label="b" fetcher={fetcher} />
      </>,
      { wrapper: wrapperFor(cache) },
    );

    await waitFor(() => expect(screen.getByTestId('a')).toHaveTextContent('rows'));
    expect(screen.getByTestId('b')).toHaveTextContent('rows');
    expect(fetcher).toHaveBeenCalledTimes(1);
  });

  it('keeps different keys apart', async () => {
    const cache = createResourceCache();
    const fetcher = vi.fn(async () => 'rows');

    render(
      <>
        <Probe label="a" fetcher={fetcher} params={{ limit: 20 }} />
        <Probe label="b" fetcher={fetcher} params={{ limit: 50 }} />
      </>,
      { wrapper: wrapperFor(cache) },
    );

    await waitFor(() => expect(fetcher).toHaveBeenCalledTimes(2));
  });

  it('does not fetch while disabled, and fetches when the gate opens', async () => {
    const cache = createResourceCache();
    const fetcher = vi.fn(async () => 'rows');

    function Gate() {
      const [enabled, setEnabled] = useState(false);
      return (
        <>
          <Probe label="a" fetcher={fetcher} enabled={enabled} />
          <button type="button" onClick={() => setEnabled(true)}>
            enable
          </button>
        </>
      );
    }

    render(<Gate />, { wrapper: wrapperFor(cache) });

    expect(screen.getByTestId('a')).toHaveTextContent('pending');
    expect(fetcher).not.toHaveBeenCalled();

    await userEvent.click(screen.getByRole('button', { name: 'enable' }));

    await waitFor(() => expect(screen.getByTestId('a')).toHaveTextContent('rows'));
    expect(fetcher).toHaveBeenCalledTimes(1);
  });

  it('refetch bypasses a fresh entry', async () => {
    const cache = createResourceCache();
    let value = 'first';
    const fetcher = vi.fn(async () => value);

    render(<Probe label="a" fetcher={fetcher} />, { wrapper: wrapperFor(cache) });
    await waitFor(() => expect(screen.getByTestId('a')).toHaveTextContent('first'));

    value = 'second';
    await userEvent.click(screen.getByRole('button', { name: 'refresh a' }));

    await waitFor(() => expect(screen.getByTestId('a')).toHaveTextContent('second'));
    expect(fetcher).toHaveBeenCalledTimes(2);
  });

  it('aborts a request whose last reader unmounted', async () => {
    const cache = createResourceCache();
    let signal: AbortSignal | undefined;
    const fetcher: ResourceFetcher<string> = ({ signal: s }) => {
      signal = s;
      return new Promise(() => {});
    };

    const view = render(<Probe label="a" fetcher={fetcher} />, { wrapper: wrapperFor(cache) });
    await waitFor(() => expect(signal).toBeDefined());
    expect(signal?.aborted).toBe(false);

    view.unmount();
    expect(signal?.aborted).toBe(true);
  });

  it('refuses to run outside a provider', () => {
    const error = vi.spyOn(console, 'error').mockImplementation(() => {});
    expect(() => render(<Probe label="a" fetcher={async () => 'rows'} />)).toThrow(
      /ResourceCacheContext/,
    );
    error.mockRestore();
  });
});
