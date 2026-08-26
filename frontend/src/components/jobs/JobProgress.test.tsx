import { screen } from '@testing-library/react';
import { describe, expect, it } from 'vitest';
import { JobProgress } from './JobProgress';
import { renderWithProviders } from '../../test/renderWithProviders';
import type { JobProgressInfo } from '../../services/contracts';

function progress(overrides: Partial<JobProgressInfo> = {}): JobProgressInfo {
  return {
    phase: 'download',
    percent: null,
    message: null,
    bytes_done: null,
    bytes_total: null,
    ...overrides,
  };
}

describe('JobProgress', () => {
  it('renders the indeterminate variant with the phase label when percent is null', () => {
    const { container } = renderWithProviders(
      <JobProgress kind="build" progress={progress({ phase: 'run', percent: null })} />,
    );
    const bar = screen.getByRole('progressbar');
    expect(bar).toHaveClass('rg-progress__bar--indeterminate');
    expect(bar).toHaveAttribute('aria-label', 'Running build');
    expect(container.querySelector('progress')).toBeNull();
  });

  it('renders a determinate bar with aria-valuenow when percent is a number', () => {
    const { container } = renderWithProviders(
      <JobProgress kind="pull" progress={progress({ percent: 38.2 })} />,
    );
    const bar = container.querySelector('progress');
    expect(bar).not.toBeNull();
    expect(bar).toHaveAttribute('aria-valuenow', '38');
  });

  it('shows a running total, and no percentage, when bytes_total is null', () => {
    renderWithProviders(
      <JobProgress
        kind="pull"
        progress={progress({ bytes_done: 1024 * 1024, bytes_total: null })}
      />,
    );
    expect(screen.getByText('1.0 MB')).toBeInTheDocument();
    expect(screen.queryByText(/%/)).not.toBeInTheDocument();
  });

  it('shows both sides when the total is known', () => {
    renderWithProviders(
      <JobProgress
        kind="pull"
        progress={progress({ bytes_done: 1024 * 1024, bytes_total: 4 * 1024 * 1024 })}
      />,
    );
    expect(screen.getByText('1.0 MB / 4.0 MB')).toBeInTheDocument();
  });

  it.each([
    ['pull', 'verify'],
    ['build', 'digest'],
    ['build', 'stage'],
  ] as const)('explains the silent phase %s/%s instead of showing a bare spinner', (kind, phase) => {
    renderWithProviders(<JobProgress kind={kind} progress={progress({ phase })} />);
    expect(screen.getByText(/silent|no output/i)).toBeInTheDocument();
  });

  it('numbers the phase so a long job reads as progress', () => {
    renderWithProviders(<JobProgress kind="pull" progress={progress({ phase: 'download' })} />);
    expect(screen.getByText(/step 4 of 8/)).toBeInTheDocument();
  });
});
