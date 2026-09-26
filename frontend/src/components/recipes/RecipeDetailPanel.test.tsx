import { describe, expect, it } from 'vitest';
import { screen } from '@testing-library/react';
import { HttpResponse, http } from 'msw';
import { RecipeDetailPanel } from './RecipeDetailPanel';
import { server } from '../../test/server';
import { API } from '../../test/handlers';
import { assetClassesFixture, recipesFixture } from '../../test/fixtures';
import { renderWithProviders } from '../../test/renderWithProviders';
import type { RecipePublic } from '../../types/api';

const FASTA_CLASS = assetClassesFixture.items.find((c) => c.name === 'fasta')!;

/**
 * A captured `fasta` recipe, given the `input_assets` blob a chained recipe
 * carries. The two recipes that output `fasta` on the public server have none,
 * and the slot name deliberately differs from the class name — that divergence
 * is the whole reason the panel renders `slot → class`.
 */
const RECIPE: RecipePublic = {
  ...recipesFixture.items[0],
  input_assets: {
    fasta_txome: { asset_class: 'fasta', description: 'Transcriptome sequences' },
  },
};

function seedClasses() {
  server.use(http.get(`${API}/asset_classes`, () => HttpResponse.json(assetClassesFixture)));
}

describe('RecipeDetailPanel', () => {
  it('renders the output asset class by name, not by id', async () => {
    seedClasses();
    renderWithProviders(<RecipeDetailPanel recipe={RECIPE} />);
    const link = await screen.findByRole('link', {
      name: `${FASTA_CLASS.name} v${FASTA_CLASS.version}`,
    });
    expect(link).toHaveAttribute('href', `/asset-classes/${RECIPE.output_asset_class_id}`);
    expect(screen.queryByText(`#${RECIPE.output_asset_class_id}`)).not.toBeInTheDocument();
  });

  it('renders each input asset as slot → linked class', async () => {
    seedClasses();
    renderWithProviders(<RecipeDetailPanel recipe={RECIPE} />);
    const link = await screen.findByRole('link', { name: FASTA_CLASS.name });
    expect(link).toHaveAttribute('href', `/asset-classes/${FASTA_CLASS.id}`);
    // The slot name is also this recipe's own name, so scope to the <code> cell.
    expect(
      screen.getAllByText('fasta_txome').some((el) => el.tagName === 'CODE'),
    ).toBe(true);
  });

  it('falls back to the raw id when the asset class index is empty', async () => {
    // The default handler returns an empty /asset_classes page.
    renderWithProviders(<RecipeDetailPanel recipe={RECIPE} />);
    const link = await screen.findByRole('link', {
      name: `#${RECIPE.output_asset_class_id}`,
    });
    expect(link).toHaveAttribute('href', `/asset-classes/${RECIPE.output_asset_class_id}`);
    // The input asset still names its class, just without a link to it.
    expect(screen.getByText('fasta')).toBeInTheDocument();
  });
});
