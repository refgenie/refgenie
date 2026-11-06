import { describe, expect, it } from 'vitest';
import { screen, waitFor, within } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { HttpResponse, http } from 'msw';
import { SpeciesPage } from './SpeciesPage';
import { server } from '../test/server';
import { API } from '../test/handlers';
import { genomesFixture } from '../test/fixtures';
import {
  localConfig,
  renderWithProviders,
  serverConfig,
} from '../test/renderWithProviders';
import type { GenomeResponse } from '../types/api';

/** A registered genome with no species recorded — it still has to be counted. */
const NAMELESS: GenomeResponse = {
  digest: '0000000000000000000000000000beef',
  aliases: ['bare-genome'],
  description: null,
  asset_count: 0,
  species_name: null,
  common_name: null,
  synonyms: null,
  taxon_id: null,
  assembly_source: null,
  assembly_accession: null,
  species: null,
};

/** A rat genome whose common name comes from the species list. */
const RAT: GenomeResponse = {
  ...NAMELESS,
  digest: '000000000000000000000000000000ab',
  aliases: ['rn7'],
  species_name: 'Rattus norvegicus',
  common_name: 'Norway rat',
  taxon_id: 10116,
  species: {
    taxon_id: 10116,
    scientific_name: 'Rattus norvegicus',
    common_name: 'Norway rat',
    synonyms: ['brown rat'],
  },
};

function seed(items: GenomeResponse[]) {
  server.use(
    http.get(`${API}/genomes`, () =>
      HttpResponse.json({
        results: items,
        pagination: { page: 0, page_size: 1000, total: items.length },
      }),
    ),
  );
}

/** The <tr> whose cells include `species`. A row's accessible name is its cells. */
async function speciesRow(species: string) {
  return screen.findByRole('row', { name: new RegExp(species) });
}

describe('SpeciesPage', () => {
  it('groups the genome list into one row per species', async () => {
    renderWithProviders(<SpeciesPage />);
    const row = await speciesRow('Homo sapiens');
    expect(within(row).getByText('human')).toBeInTheDocument();
    expect(within(row).getByText('9606')).toBeInTheDocument();
    // The two fixture genomes carry 3 and 1 assets.
    expect(within(row).getByText('2')).toBeInTheDocument();
    expect(within(row).getByText('4')).toBeInTheDocument();
  });

  it('counts the species and the ones with assets in the lede', async () => {
    renderWithProviders(<SpeciesPage />);
    expect(
      await screen.findByText(/1 species · 1 with assets built/),
    ).toBeInTheDocument();
  });

  it('links a row into the genome list filtered to that species', async () => {
    renderWithProviders(<SpeciesPage />);
    expect(await screen.findByRole('link', { name: 'Homo sapiens' })).toHaveAttribute(
      'href',
      '/genomes?q=Homo+sapiens&fields=species&op=contains',
    );
  });

  it("matches the species list's other names without showing them", async () => {
    seed([...genomesFixture.results, RAT]);
    renderWithProviders(<SpeciesPage />);
    await screen.findByText('Homo sapiens');
    expect(screen.queryByText(/brown rat/)).not.toBeInTheDocument();
    await userEvent.type(screen.getByLabelText('Search species'), 'brown rat');
    const row = await speciesRow('Rattus norvegicus');
    expect(within(row).getByText('Norway rat')).toBeInTheDocument();
    await waitFor(() =>
      expect(screen.queryByText('Homo sapiens')).not.toBeInTheDocument(),
    );
    expect(screen.getByRole('table')).not.toHaveTextContent('brown rat');
  });

  it('shows the search empty state when nothing matches', async () => {
    renderWithProviders(<SpeciesPage />);
    await screen.findByText('Homo sapiens');
    await userEvent.type(screen.getByLabelText('Search species'), 'mouse');
    expect(await screen.findByText(/No results match/)).toBeInTheDocument();
  });

  it('renders a genome with no species as an unlinked Unspecified row', async () => {
    seed([...genomesFixture.results, NAMELESS]);
    renderWithProviders(<SpeciesPage />);
    expect(await screen.findByText('Unspecified')).toBeInTheDocument();
    expect(screen.queryByRole('link', { name: 'Unspecified' })).not.toBeInTheDocument();
  });

  it('offers the tree link on a server, where the tree is reachable', async () => {
    renderWithProviders(<SpeciesPage />, { config: serverConfig });
    expect(await screen.findByRole('link', { name: 'Show in tree' })).toHaveAttribute(
      'href',
      '/tree?species=Homo%20sapiens',
    );
  });

  it('has no tree column on a local dash, where /tree is not offered', async () => {
    renderWithProviders(<SpeciesPage />, { config: localConfig });
    await screen.findByText('Homo sapiens');
    expect(
      screen.queryByRole('link', { name: 'Show in tree' }),
    ).not.toBeInTheDocument();
  });
});
