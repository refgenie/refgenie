import { screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { HttpResponse, http } from 'msw';
import { beforeEach, describe, expect, it } from 'vitest';
import { PullButton } from './PullButton';
import { useJobStore } from '../../stores/jobStore';
import { renderWithProviders } from '../../test/renderWithProviders';
import { server } from '../../test/server';
import { LOCAL_API } from '../../test/handlers';
import { makeJob } from '../../test/jobFixtures';

const props = {
  serverUrl: 'https://refgenomes.databio.org',
  genomeDigest: 'genome-digest-1',
  assetGroupName: 'fasta',
  assetName: 'default',
  assetDigest: 'asset-digest-1',
  archiveDigest: 'archive-1',
  archiveSize: 1024,
  existsLocally: false,
};

beforeEach(() => useJobStore.getState().reset());

describe('PullButton', () => {
  it('replaces the button with inline progress when a job is already running', () => {
    useJobStore.getState().upsertJob(
      makeJob({
        id: 'j-live',
        status: 'running',
        progress: {
          phase: 'download',
          percent: 12,
          message: null,
          bytes_done: null,
          bytes_total: null,
        },
      }),
    );

    renderWithProviders(<PullButton {...props} />);

    expect(screen.queryByRole('button', { name: 'Pull' })).not.toBeInTheDocument();
    expect(screen.getByRole('button', { name: 'Show console' })).toBeInTheDocument();
    expect(screen.getByText('Downloading')).toBeInTheDocument();
  });

  it('sends the action header and registers the returned job', async () => {
    let headerSeen: string | null = null;
    let body: unknown;
    server.use(
      http.post(`${LOCAL_API}/actions/pull`, async ({ request }) => {
        headerSeen = request.headers.get('x-refgenie-action');
        body = await request.json();
        return HttpResponse.json(
          {
            job_id: 'j-new',
            kind: 'pull',
            status: 'queued',
            created_at: '2026-08-12T10:00:00Z',
            duplicate: false,
          },
          { status: 202 },
        );
      }),
    );

    renderWithProviders(<PullButton {...props} />);
    await userEvent.click(screen.getByRole('button', { name: 'Pull' }));

    await waitFor(() => expect(useJobStore.getState().jobs['j-new']).toBeDefined());
    expect(headerSeen).not.toBeNull();
    // Never omitted: the server default exists so the puller cannot reach a
    // stdin prompt, and a UI relying on it is one refactor from hanging a worker.
    expect(body).toMatchObject({
      force: false,
      asset_group: 'fasta',
      asset: 'default',
      genome_digest: 'genome-digest-1',
    });
    // Exactly one genome reference: sending both is a 422.
    expect(body).not.toHaveProperty('genome');
  });

  it('focuses the existing card on a duplicate instead of creating a second one', async () => {
    server.use(
      http.post(`${LOCAL_API}/actions/pull`, () =>
        HttpResponse.json(
          {
            job_id: 'j-existing',
            kind: 'pull',
            status: 'running',
            created_at: '2026-08-12T10:00:00Z',
            duplicate: true,
          },
          { status: 202 },
        ),
      ),
    );

    renderWithProviders(<PullButton {...props} />);
    await userEvent.click(screen.getByRole('button', { name: 'Pull' }));

    await waitFor(() =>
      expect(useJobStore.getState().focusedJobId).toBe('j-existing'),
    );
    // No second card, and no error surface: a duplicate is a 202, not a 409.
    expect(Object.keys(useJobStore.getState().jobs)).toHaveLength(0);
    expect(screen.queryByRole('alert')).not.toBeInTheDocument();
  });

  it('disables the primary button and offers an overwrite when the asset is local', async () => {
    renderWithProviders(<PullButton {...props} existsLocally />);
    expect(screen.getByRole('button', { name: 'Already local' })).toBeDisabled();

    await userEvent.click(screen.getByRole('button', { name: 'More pull options' }));
    expect(screen.getByRole('menuitem', { name: 'Re-pull (overwrite)' })).toBeInTheDocument();
  });
});
