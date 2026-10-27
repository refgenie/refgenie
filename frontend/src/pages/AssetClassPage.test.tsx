import { describe, expect, it } from 'vitest';
import { screen } from '@testing-library/react';
import { HttpResponse, http } from 'msw';
import { Route, Routes } from 'react-router-dom';
import { AssetClassPage } from './AssetClassPage';
import { server } from '../test/server';
import { API } from '../test/handlers';
import { assetClassesFixture, assetGroupsFixture, recipesFixture } from '../test/fixtures';
import { renderWithProviders } from '../test/renderWithProviders';
import type { AssetGroupPublic, GenomeResponse } from '../types/api';
import type { Paginated } from '../types/pagination';

/** `fasta` on the public server; the fixtures were captured for it. */
const FASTA_CLASS = assetClassesFixture.items.find((c) => c.name === 'fasta')!;
const GROUPS = assetGroupsFixture.items;

function genome(digest: string, alias: string, species: string): GenomeResponse {
  return {
    digest,
    aliases: [alias],
    description: null,
    asset_count: 1,
    species_name: species,
    common_name: null,
    taxon_id: null,
    assembly_source: null,
    assembly_accession: null,
  };
}

/** One named genome per captured asset group, so every row can resolve. */
const GENOME_INDEX: Paginated<GenomeResponse> = {
  items: GROUPS.map((group, i) => genome(group.genome_digest, `genome${i}`, 'Homo sapiens')),
  pagination: { offset: 0, limit: 1000, total: GROUPS.length },
};

function page<T>(items: T[]): Paginated<T> {
  return { items, pagination: { offset: 0, limit: 50, total: items.length } };
}

function seed(groups: AssetGroupPublic[] = GROUPS) {
  server.use(
    http.get(`${API}/asset_classes/:id`, () => HttpResponse.json(FASTA_CLASS)),
    http.get(`${API}/asset_classes`, () => HttpResponse.json(assetClassesFixture)),
    http.get(`${API}/asset_groups`, () => HttpResponse.json(page(groups))),
    http.get(`${API}/genomes`, () => HttpResponse.json(GENOME_INDEX)),
    http.get(`${API}/recipes`, () => HttpResponse.json(recipesFixture)),
  );
}

function renderClass() {
  return renderWithProviders(
    <Routes>
      <Route path="/asset-classes/:id" element={<AssetClassPage />} />
    </Routes>,
    { route: `/asset-classes/${FASTA_CLASS.id}` },
  );
}

describe('AssetClassPage', () => {
  it('renders a genome alias, not a bare digest, in the Genome column', async () => {
    seed();
    renderClass();
    const link = await screen.findByRole('link', { name: 'genome0' });
    expect(link).toHaveAttribute('href', `/genomes/${GROUPS[0].genome_digest}`);
  });

  it('deep-links each row to its asset group', async () => {
    seed();
    renderClass();
    await screen.findByRole('link', { name: 'genome0' });
    const groupLinks = screen.getAllByRole('link', { name: GROUPS[0].name });
    expect(groupLinks.map((a) => a.getAttribute('href'))).toContain(
      `/asset-groups/${GROUPS[0].id}`,
    );
  });

  it('lists the recipes that build the class under "Built by"', async () => {
    seed();
    renderClass();
    const recipe = recipesFixture.items[0];
    const link = await screen.findByRole('link', { name: recipe.name });
    expect(link).toHaveAttribute('href', `/recipes/${recipe.id}`);
  });

  it('drops an asset group belonging to a sibling version of the class', async () => {
    const sibling: AssetGroupPublic = {
      ...GROUPS[0],
      id: 9999,
      asset_class_id: FASTA_CLASS.id! + 1000,
      genome_digest: GROUPS[0].genome_digest,
    };
    seed([sibling]);
    renderClass();
    expect(
      await screen.findByText('No genome on this server has this asset class built.'),
    ).toBeInTheDocument();
  });
});
