import { describe, expect, it, vi } from 'vitest';
import { act, render, screen } from '@testing-library/react';
import { useEffect, useState } from 'react';
import { createMemoryRouter, RouterProvider, useParams } from 'react-router-dom';

const unmounted: string[] = [];

vi.mock('../pages/GenomePage', () => ({
  // State taken on mount, as the page's download, extract and compare state is.
  GenomePage: function FakeGenomePage() {
    const { digest } = useParams<{ digest: string }>();
    const [first] = useState(digest);
    useEffect(() => () => void unmounted.push(first ?? ''), [first]);
    return <p>state for {first}</p>;
  },
}));

const { GenomeRoute } = await import('./routes');

describe('GenomeRoute', () => {
  it('remounts the page when the digest changes, dropping its state', async () => {
    const router = createMemoryRouter(
      [{ path: '/genomes/:digest', element: <GenomeRoute /> }],
      { initialEntries: ['/genomes/AAA'] },
    );
    render(<RouterProvider router={router} />);
    expect(screen.getByText('state for AAA')).toBeInTheDocument();
    await act(() => router.navigate('/genomes/BBB?tab=download'));
    expect(screen.getByText('state for BBB')).toBeInTheDocument();
    // The old page unmounted, so its cleanups (aborting downloads) ran.
    expect(unmounted).toEqual(['AAA']);
  });
});
