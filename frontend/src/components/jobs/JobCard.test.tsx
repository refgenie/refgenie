import { screen } from '@testing-library/react';
import { beforeEach, describe, expect, it } from 'vitest';
import { JobCard } from './JobCard';
import { useJobStore } from '../../stores/jobStore';
import { renderWithProviders } from '../../test/renderWithProviders';
import { makeJob } from '../../test/jobFixtures';

beforeEach(() => useJobStore.getState().reset());

describe('JobCard', () => {
  it('shows the queue position and nothing else for a queued job', () => {
    const job = makeJob({ status: 'queued', queue_position: 2 });
    const { container } = renderWithProviders(<JobCard job={job} />);

    expect(screen.getByText('Queued — 2 ahead')).toBeInTheDocument();
    // A queued card must not show a stalled-looking bar or an empty log box.
    expect(container.querySelector('.rg-progress')).toBeNull();
    expect(container.querySelector('.rg-log-tail')).toBeNull();
  });

  it('says "next up" rather than "0 ahead"', () => {
    renderWithProviders(<JobCard job={makeJob({ status: 'queued', queue_position: 0 })} />);
    expect(screen.getByText('Queued — next up')).toBeInTheDocument();
  });

  it('says it is waiting for output when a running job has no log lines yet', () => {
    renderWithProviders(<JobCard job={makeJob({ status: 'running' })} />);
    expect(screen.getByText(/waiting for output/i)).toBeInTheDocument();
  });

  it('renders a build failure as a preformatted log with a View log action', () => {
    const job = makeJob({
      kind: 'build',
      status: 'failed',
      error: {
        code: 'build_failed',
        message: 'Build failed',
        detail: 'Command exited with 1\nbowtie2-build: not found',
        field: null,
      },
    });
    const { container } = renderWithProviders(<JobCard job={job} />);

    expect(screen.getByText('Build failed')).toBeInTheDocument();
    const block = container.querySelector('pre.rg-job-error__log');
    expect(block?.textContent).toContain('bowtie2-build: not found');
    expect(screen.getByRole('button', { name: 'View log' })).toBeInTheDocument();
  });

  it('disables cancel for a running build rather than lying about it', () => {
    const job = makeJob({ kind: 'build', status: 'running', cancellable: false });
    renderWithProviders(<JobCard job={job} />);
    expect(screen.getByRole('button', { name: 'Cancel' })).toBeDisabled();
  });

  it('says nothing was rebuilt when the asset already existed', () => {
    const job = makeJob({ kind: 'build', status: 'succeeded', skipped: true });
    renderWithProviders(<JobCard job={job} />);
    expect(screen.getByText(/already built/i)).toBeInTheDocument();
  });
});
