import { describe, expect, it } from 'vitest';
import { emptyBuildForm, registryPathPreview, validateBuildForm } from './validate';
import type { RecipePublic } from '../../types/api';

const recipe: RecipePublic = {
  id: 1,
  name: 'bowtie2_index',
  version: '0.0.1',
  description: 'Bowtie2 index',
  output_asset_class_id: 2,
  command_templates: ['bowtie2-build ...'],
  // One required file, one defaulted param, one required input asset: the three
  // shapes the form has to generate from.
  input_files: { fasta: { description: 'A FASTA file' } },
  input_params: { threads: { description: 'CPU threads', default: 4 } },
  input_assets: { parent: { asset_class: 'fasta', description: 'The fasta asset' } },
  docker_image: null,
  custom_seek_keys: null,
  default_asset: '{genome}',
  inherent: null,
};

function values(overrides: Partial<ReturnType<typeof emptyBuildForm>> = {}) {
  return {
    ...emptyBuildForm(),
    genome: 'hg38',
    recipeName: 'bowtie2_index',
    recipeVersion: '0.0.1',
    assetGroupName: 'bowtie2_index',
    ...overrides,
  };
}

describe('validateBuildForm', () => {
  it('requires a genome, a recipe and an asset group', () => {
    const errors = validateBuildForm(null, emptyBuildForm());
    expect(errors.genome).toBeTruthy();
    expect(errors.recipe).toBeTruthy();
    expect(errors.asset_group).toBeTruthy();
  });

  it('flags a required file and a required input asset, but not a defaulted param', () => {
    const errors = validateBuildForm(recipe, values());
    expect(errors['params.files.fasta']).toBeTruthy();
    expect(errors['params.assets.parent']).toBeTruthy();
    expect(errors['params.params.threads']).toBeUndefined();
  });

  it('passes once every required field is filled', () => {
    const errors = validateBuildForm(
      recipe,
      values({
        files: { fasta: '/data/hg38.fa' },
        assets: { parent: 'hg38/fasta:default' },
      }),
    );
    expect(errors).toEqual({});
  });

  it('rejects registry-path characters the backend blacklists', () => {
    for (const bad of ['my group', 'my:group', 'my/group']) {
      const errors = validateBuildForm(recipe, values({ assetGroupName: bad }));
      expect(errors.asset_group, bad).toBeTruthy();
    }
  });
});

describe('registryPathPreview', () => {
  it('says "auto" when the asset name is left to the recipe', () => {
    expect(registryPathPreview(values())).toBe('hg38/bowtie2_index:auto');
  });
});
