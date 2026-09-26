import { describe, expect, it } from 'vitest';
import { render, screen } from '@testing-library/react';
import { ErrorState } from './states';
import { ApiError } from '../../services/http';

const notFound = (code: string, detail: string) =>
  new ApiError({ status: 404, detail, url: '/v4/genomes/abc', code });

describe('ErrorState on a 404', () => {
  it('shows only the identifier when the server had nothing to add', () => {
    render(
      <ErrorState error={notFound('not_found', 'Genome abc not found')} subject="abc" />,
    );
    expect(screen.getByText('Not found:')).toBeInTheDocument();
    expect(screen.queryByText(/Genome abc not found/)).not.toBeInTheDocument();
  });

  it('shows the server message when the code carries an explanation', () => {
    // Swallowing this is what left a stale alias looking like a plain 404 with
    // no hint of which half of the instance was wrong.
    const message =
      'No genome record for abc, but the alias(es) rCRSd still point at it. ' +
      'Remove the stale alias with `refgenie alias remove rCRSd`.';
    render(<ErrorState error={notFound('stale_alias', message)} subject="abc" />);
    expect(screen.getByText(message)).toBeInTheDocument();
    // The copyable identifier survives alongside it.
    expect(screen.getByText('Not found:')).toBeInTheDocument();
  });
});
