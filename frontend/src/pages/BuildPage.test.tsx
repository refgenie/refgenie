import { screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { HttpResponse, http } from 'msw';
import { beforeEach, describe, expect, it } from 'vitest';
import { BuildPage } from './BuildPage';
import { useJobStore } from '../stores/jobStore';
import { renderWithProviders } from '../test/renderWithProviders';
import { server } from '../test/server';
import { API, LOCAL_API } from '../test/handlers';
import type { RecipePublic } from '../types/api';

const recipe: RecipePublic = {
  id: 1,
  name: 'bowtie2_index',
  version: '0.0.1',
  description: 'Bowtie2 index',
  output_asset_class_id: 2,
  command_templates: [],
  input_params: {},
  input_files: { fasta: { description: 'A FASTA file' } },
  input_assets: {},
  docker_image: null,
  custom_seek_keys: null,
  default_asset: '{genome}',
  inherent: null,
};

function seedCatalog() {
  server.use(
    http.get(`${API}/recipes`, () =>
      HttpResponse.json({
        items: [recipe],
        pagination: { offset: 0, limit: 200, total: 1 },
      }),
    ),
    http.get(`${API}/asset_classes`, () =>
      HttpResponse.json({
        items: [
          { id: 2, name: 'bowtie2_index', version: '0.0.1', description: null, serving_modes: [] },
        ],
        pagination: { offset: 0, limit: 200, total: 1 },
      }),
    ),
    http.get(`${API}/genomes`, () =>
      HttpResponse.json({
        items: [
          {
            digest: 'genome-digest-1',
            aliases: ['hg38'],
            description: null,
            asset_count: 1,
            species_name: null,
            common_name: null,
            taxon_id: null,
            assembly_source: null,
            assembly_accession: null,
          },
        ],
        pagination: { offset: 0, limit: 25, total: 1 },
      }),
    ),
    http.get(`${API}/configurations`, () =>
      HttpResponse.json({
        items: [{ version: 1, servers: [], genome_folder: '/data', genome_stage_folder: null }],
        pagination: { offset: 0, limit: 1, total: 1 },
      }),
    ),
  );
}

beforeEach(() => useJobStore.getState().reset());

describe('BuildPage', () => {
  it('submits the expected body with the action header and registers the job', async () => {
    seedCatalog();
    let headerSeen: string | null = null;
    let body: Record<string, unknown> | undefined;
    server.use(
      http.post(`${LOCAL_API}/actions/build/preflight`, () =>
        HttpResponse.json({ ok: true, errors: [], resolved: { asset_name: 'hg38' } }),
      ),
      http.post(`${LOCAL_API}/actions/build`, async ({ request }) => {
        headerSeen = request.headers.get('x-refgenie-action');
        body = (await request.json()) as Record<string, unknown>;
        return HttpResponse.json(
          {
            job_id: 'j-build',
            kind: 'build',
            status: 'queued',
            created_at: '2026-08-12T10:00:00Z',
            duplicate: false,
          },
          { status: 202 },
        );
      }),
    );

    renderWithProviders(<BuildPage />, { route: '/build' });

    await userEvent.type(await screen.findByLabelText(/^Genome/), 'hg38');
    await userEvent.selectOptions(
      await screen.findByLabelText(/^Recipe \*/),
      'bowtie2_index',
    );
    await userEvent.type(screen.getByLabelText(/fasta/), '/data/hg38.fa');

    await userEvent.click(screen.getByRole('button', { name: 'Build' }));

    await waitFor(() => expect(useJobStore.getState().jobs['j-build']).toBeDefined());
    expect(headerSeen).not.toBeNull();
    // Backend field names, with the recipe inputs NESTED under `params`; the
    // models are extra="forbid", so a stale name here would be a 422.
    expect(body).toMatchObject({
      genome: 'hg38',
      asset_group: 'bowtie2_index',
      recipe: 'bowtie2_index',
      params: { files: { fasta: '/data/hg38.fa' } },
    });
    expect(body).not.toHaveProperty('docker');
    expect(body).not.toHaveProperty('asset_group_name');
  });

  it('renders a preflight field error under the input it names', async () => {
    seedCatalog();
    server.use(
      http.post(`${LOCAL_API}/actions/build/preflight`, () =>
        HttpResponse.json({
          ok: false,
          errors: [
            {
              field: 'params.files.fasta',
              code: 'missing_build_input',
              message: 'No such file',
            },
          ],
          resolved: {},
        }),
      ),
    );

    renderWithProviders(<BuildPage />, { route: '/build' });

    await userEvent.type(await screen.findByLabelText(/^Genome/), 'hg38');
    await userEvent.selectOptions(
      await screen.findByLabelText(/^Recipe \*/),
      'bowtie2_index',
    );
    await userEvent.type(screen.getByLabelText(/fasta/), '/nope.fa');

    expect(await screen.findByText('No such file')).toBeInTheDocument();
  });

  it('surfaces a coarse preflight error that no single input owns', async () => {
    seedCatalog();
    server.use(
      http.post(`${LOCAL_API}/actions/build/preflight`, () =>
        HttpResponse.json({
          ok: false,
          errors: [
            { field: 'params.assets', code: 'asset_not_found', message: 'No fasta asset' },
          ],
          resolved: {},
        }),
      ),
    );

    renderWithProviders(<BuildPage />, { route: '/build' });
    await userEvent.type(await screen.findByLabelText(/^Genome/), 'hg38');
    await userEvent.selectOptions(await screen.findByLabelText(/^Recipe \*/), 'bowtie2_index');
    await userEvent.type(screen.getByLabelText(/fasta/), '/data/hg38.fa');

    expect(await screen.findByText('No fasta asset')).toBeInTheDocument();
  });

  it('hides the staging checkbox when the configuration has no stage folder', async () => {
    seedCatalog();
    renderWithProviders(<BuildPage />, { route: '/build' });
    await screen.findByLabelText(/^Recipe \*/);
    // Offering it unconditionally manufactures a bare ValueError server-side.
    expect(screen.queryByLabelText(/Stage the asset/)).not.toBeInTheDocument();
  });
});
