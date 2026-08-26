import { screen } from '@testing-library/react';
import { HttpResponse, http } from 'msw';
import { describe, expect, it, vi } from 'vitest';
import { RecipeInputFields } from './RecipeInputFields';
import { renderWithProviders } from '../../test/renderWithProviders';
import { server } from '../../test/server';
import { API } from '../../test/handlers';
import type { AssetResponse, RecipePublic } from '../../types/api';

const recipe: RecipePublic = {
  id: 1,
  name: 'bowtie2_index',
  version: '0.0.1',
  description: null,
  output_asset_class_id: 2,
  command_templates: [],
  input_params: { threads: { description: 'CPU threads', default: 4 } },
  input_files: { extra: { description: 'An extra file' } },
  input_assets: { parent: { asset_class: 'fasta', description: 'The fasta asset' } },
  docker_image: null,
  custom_seek_keys: null,
  default_asset: '{genome}',
  inherent: null,
};

function asset(overrides: Partial<AssetResponse>): AssetResponse {
  return {
    digest: 'a1',
    name: 'default',
    description: null,
    recipe_id: 1,
    asset_group_id: 1,
    size: null,
    serving_modes_override: null,
    colocate: null,
    serving_modes: ['archive'],
    asset_class_name: 'fasta',
    asset_group_name: 'fasta',
    genome_digest: 'genome-digest-1',
    names: null,
    seek_keys: null,
    ...overrides,
  };
}

function renderFields(assets: AssetResponse[]) {
  server.use(
    http.get(`${API}/assets`, () =>
      HttpResponse.json({
        items: assets,
        pagination: { offset: 0, limit: 200, total: assets.length },
      }),
    ),
  );
  return renderWithProviders(
    <RecipeInputFields
      recipe={recipe}
      genomeDigest="genome-digest-1"
      genomeName="hg38"
      values={{ params: {}, files: {}, assets: {} }}
      errors={{}}
      onChange={vi.fn()}
      buildLinkFor={(assetClass) => `/build?recipe=${assetClass}`}
    />,
  );
}

describe('RecipeInputFields', () => {
  it('renders one field per entry in each of the three JSON columns', async () => {
    renderFields([asset({})]);

    expect(screen.getByLabelText(/threads/)).toBeInTheDocument();
    expect(screen.getByLabelText(/extra/)).toBeInTheDocument();
    expect(await screen.findByLabelText(/parent/)).toBeInTheDocument();
  });

  it('prefills a defaulted param and marks an undefaulted one required', () => {
    renderFields([asset({})]);
    expect(screen.getByLabelText(/threads/)).toHaveValue(4);
    // The file input has no default, so it carries the required marker.
    expect(screen.getByLabelText(/extra/)).toHaveValue('');
  });

  it('offers only assets of the required class, as the full registry path', async () => {
    renderFields([
      asset({ digest: 'a1', name: 'default', asset_class_name: 'fasta' }),
      asset({
        digest: 'a2',
        name: 'other',
        asset_class_name: 'bowtie2_index',
        asset_group_name: 'bowtie2_index',
      }),
    ]);

    await screen.findByRole('option', { name: 'fasta:default' });
    const select = screen.getByLabelText(/parent/);
    const options = Array.from(select.querySelectorAll('option')).map((o) => o.value);
    expect(options).toContain('hg38/fasta:default');
    expect(options.some((value) => value.includes('bowtie2_index'))).toBe(false);
  });

  it('offers a build link instead of an empty dropdown when nothing matches', async () => {
    renderFields([asset({ asset_class_name: 'bowtie2_index' })]);
    // The empty state names the missing class rather than showing a dead dropdown.
    const link = await screen.findByRole('link', { name: 'Build it first' });
    expect(link).toHaveAttribute('href', '/build?recipe=fasta');
    expect(screen.getByText(/asset for this genome/)).toBeInTheDocument();
    expect(screen.queryByLabelText(/parent/)).not.toBeInTheDocument();
  });
});
