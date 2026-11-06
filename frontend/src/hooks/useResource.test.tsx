import { StrictMode } from 'react';
import { describe, expect, it, vi } from 'vitest';
import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { ResourceCacheContext, useResource } from './useResource';
import { createResourceCache } from '../services/resourceCache';

describe('useResource', () => {
  it('renders data under StrictMode and refetches on demand', async () => {
    let value = 'first';
    const fetcher = vi.fn(async () => value);

    function Probe() {
      const query = useResource(['genomes', {}], fetcher, { staleTime: 60_000 });
      return (
        <button type="button" onClick={() => void query.refetch()}>
          {query.isPending ? 'pending' : query.data}
        </button>
      );
    }

    render(
      <StrictMode>
        <ResourceCacheContext.Provider value={createResourceCache()}>
          <Probe />
        </ResourceCacheContext.Provider>
      </StrictMode>,
    );

    const button = await screen.findByRole('button', { name: 'first' });
    value = 'second';
    await userEvent.click(button);
    expect(await screen.findByRole('button', { name: 'second' })).toBeInTheDocument();
  });
});
