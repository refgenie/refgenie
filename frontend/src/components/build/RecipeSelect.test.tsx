import { screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { HttpResponse, http } from 'msw';
import { describe, expect, it, vi } from 'vitest';
import { RecipeSelect } from './RecipeSelect';
import { renderWithProviders } from '../../test/renderWithProviders';
import { server } from '../../test/server';
import { API } from '../../test/handlers';
import type { RecipePublic } from '../../types/api';

function recipe(overrides: Partial<RecipePublic> & { name: string; id: number }): RecipePublic {
  return {
    version: '0.0.1',
    description: null,
    output_asset_class_id: 1,
    command_templates: [],
    input_params: {},
    input_files: {},
    input_assets: {},
    docker_image: null,
    custom_seek_keys: null,
    default_asset: '{genome}',
    inherent: null,
    ...overrides,
  };
}

/*
 * Mirrors the shape of the installed catalog that made this list read as
 * doubled: one class with two recipes (`fasta`), one class whose single recipe
 * is named differently (`bed` -> `bed12`), and one where class and recipe agree.
 */
const RECIPES: RecipePublic[] = [
  recipe({ id: 1, name: 'fasta', output_asset_class_id: 10 }),
  recipe({ id: 2, name: 'fasta_txome', output_asset_class_id: 10 }),
  recipe({ id: 3, name: 'bed12', output_asset_class_id: 11 }),
  recipe({ id: 4, name: 'bowtie2_index', output_asset_class_id: 12, version: '0.0.1' }),
  recipe({ id: 5, name: 'bowtie2_index', output_asset_class_id: 12, version: '0.0.2' }),
];

const CLASSES = [
  { id: 10, name: 'fasta', version: '0.0.1', description: null, serving_modes: [] },
  { id: 11, name: 'bed', version: '0.0.1', description: null, serving_modes: [] },
  { id: 12, name: 'bowtie2_index', version: '0.0.1', description: null, serving_modes: [] },
];

/** Distinct recipe names: versions collapse into one row with a version field. */
const RECIPE_COUNT = new Set(RECIPES.map((item) => item.name)).size;

function seedCatalog() {
  server.use(
    http.get(`${API}/recipes`, () =>
      HttpResponse.json({
        items: RECIPES,
        pagination: { offset: 0, limit: 200, total: RECIPES.length },
      }),
    ),
    http.get(`${API}/asset_classes`, () =>
      HttpResponse.json({
        items: CLASSES,
        pagination: { offset: 0, limit: 200, total: CLASSES.length },
      }),
    ),
  );
}

function renderSelect(props: Partial<Parameters<typeof RecipeSelect>[0]> = {}) {
  seedCatalog();
  return renderWithProviders(
    <RecipeSelect
      recipeName={null}
      recipeVersion={null}
      onChange={vi.fn()}
      {...props}
    />,
  );
}

async function recipeSelect() {
  const field = await screen.findByLabelText(/^Recipe \*/);
  await waitFor(() => expect(field.querySelectorAll('option').length).toBeGreaterThan(1));
  return field as HTMLSelectElement;
}

describe('RecipeSelect', () => {
  it('renders exactly one selectable row per recipe and no group headers', async () => {
    renderSelect();
    const field = await recipeSelect();

    // The invariant: every row is a recipe the user can choose.
    expect(field.querySelectorAll('optgroup')).toHaveLength(0);
    expect([...field.children].every((child) => child.tagName === 'OPTION')).toBe(true);

    const options = [...field.querySelectorAll('option')];
    expect(options).toHaveLength(RECIPE_COUNT + 1); // + the placeholder
    expect(options.filter((option) => option.disabled)).toHaveLength(0);
  });

  it('names the asset class inline only where it differs from the recipe', async () => {
    renderSelect();
    const field = await recipeSelect();

    expect([...field.querySelectorAll('option')].map((option) => option.textContent)).toEqual([
      'Choose a recipe…',
      'bed · bed12',
      'bowtie2_index',
      'fasta',
      'fasta · fasta_txome',
    ]);
  });

  it('reports the picked recipe by name, not by its label', async () => {
    const onChange = vi.fn();
    renderSelect({ onChange });
    const field = await recipeSelect();

    await userEvent.selectOptions(field, 'fasta_txome');

    expect(onChange).toHaveBeenCalledWith(expect.objectContaining({ name: 'fasta_txome' }));
  });

  it('still offers the version sub-field for a recipe with several versions', async () => {
    renderSelect({ recipeName: 'bowtie2_index', recipeVersion: '0.0.2' });
    await recipeSelect();

    const versions = (await screen.findByLabelText(/^Recipe version/)) as HTMLSelectElement;
    expect([...versions.options].map((option) => option.value)).toEqual(['0.0.2', '0.0.1']);
    expect(versions.value).toBe('0.0.2');
  });
});
