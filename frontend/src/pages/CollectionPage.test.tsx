import { beforeAll, describe, expect, it } from 'vitest';
import { screen, waitFor, within } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { HttpResponse, http } from 'msw';
import { Route, Routes } from 'react-router-dom';
import { CollectionPage } from './CollectionPage';
import { server } from '../test/server';
import { API } from '../test/handlers';
import { renderWithProviders, serverConfig } from '../test/renderWithProviders';
import { vgpCollection, vgpDetail, vgpGenomes } from '../test/collectionFixtures';

beforeAll(() => {
  // jsdom does not implement ResizeObserver, and the tree observes its panel.
  globalThis.ResizeObserver = class {
    observe() {}
    unobserve() {}
    disconnect() {}
  } as unknown as typeof ResizeObserver;
});

const config = { ...serverConfig, collections: [vgpCollection] };

function seed() {
  server.use(
    http.get(`${API}/genomes`, ({ request }) => {
      const store = new URL(request.url).searchParams.get('store');
      if (store !== 'vgp')
        return HttpResponse.json({ detail: 'bad store' }, { status: 422 });
      return HttpResponse.json({
        results: vgpGenomes,
        pagination: { page: 0, page_size: 1000, total: vgpGenomes.length },
      });
    }),
    http.get(`${API}/genomes/:digest`, () => HttpResponse.json(vgpDetail)),
  );
}

function renderAt(route: string) {
  return renderWithProviders(
    <Routes>
      <Route path="/collections" element={<p>collections index</p>} />
      <Route path="/collections/:id" element={<CollectionPage />} />
    </Routes>,
    { config, route },
  );
}

async function loaded() {
  await screen.findByRole('link', { name: 'fLatCha1.pri' });
}

function rows() {
  return screen.getAllByRole('row').slice(1);
}

describe('CollectionPage', () => {
  it('shows the branding: name, blurb, links and note', async () => {
    seed();
    renderAt('/collections/vgp');
    expect(
      screen.getByRole('heading', { level: 1, name: vgpCollection.name }),
    ).toBeInTheDocument();
    expect(screen.getByText(vgpCollection.blurb)).toBeInTheDocument();
    for (const link of vgpCollection.links) {
      expect(
        screen.getByRole('link', { name: new RegExp(link.label) }),
      ).toHaveAttribute('href', link.url);
    }
    expect(screen.getByText(vgpCollection.note!)).toBeInTheDocument();
    await loaded();
  });

  it('shows the placeholder badge and no logo image', async () => {
    seed();
    const { container } = renderAt('/collections/vgp');
    expect(screen.getByText('VGP')).toBeInTheDocument();
    expect(container.querySelector('img.rg-collection-logo')).toBeNull();
    await loaded();
  });

  it('shows headline counts', async () => {
    seed();
    renderAt('/collections/vgp');
    await loaded();
    const value = (term: string) =>
      screen.getByText(term, { selector: 'dt' }).nextElementSibling?.textContent;
    expect(value('Genomes')).toBe('3');
    expect(value('Species')).toBe('2');
    expect(value('Genomes with built assets')).toBe('1');
  });

  it('shows accession and haplotype columns', async () => {
    seed();
    renderAt('/collections/vgp');
    await loaded();
    const accession = screen.getByRole('link', { name: /GCF_902713425\.1/ });
    expect(accession).toHaveAttribute(
      'href',
      'https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_902713425.1/',
    );
    const row = accession.closest('tr')!;
    expect(within(row).getByText('maternal')).toBeInTheDocument();
    expect(within(row).getByText('Chromosome')).toBeInTheDocument();
  });

  it('filters by search text and by haplotype', async () => {
    seed();
    renderAt('/collections/vgp');
    await loaded();
    expect(rows()).toHaveLength(3);
    await userEvent.type(
      screen.getByLabelText('Search this collection'),
      'GCA_905165185',
    );
    await waitFor(() => expect(rows()).toHaveLength(1));
    expect(
      screen.getByRole('link', { name: 'aBomVar4.paternal.1' }),
    ).toBeInTheDocument();

    await userEvent.clear(screen.getByLabelText('Search this collection'));
    await userEvent.selectOptions(screen.getByLabelText('Haplotype'), 'maternal');
    await waitFor(() => expect(rows()).toHaveLength(1));
    expect(
      screen.getByRole('link', { name: 'fAciRut3.2_maternal_haplotype' }),
    ).toBeInTheDocument();
  });

  it('opens a detail panel that says when no quality metrics exist', async () => {
    seed();
    renderAt('/collections/vgp');
    await loaded();
    await userEvent.click(
      screen.getByRole('button', { name: 'Details for fAciRut3.2_maternal_haplotype' }),
    );
    const dialog = await screen.findByRole('heading', {
      level: 2,
      name: 'fAciRut3.2_maternal_haplotype',
    });
    expect(dialog).toBeInTheDocument();
    expect(
      await screen.findByText(
        'No quality metrics have been published for this assembly.',
      ),
    ).toBeInTheDocument();
    const modal = document.querySelector<HTMLElement>('.modal')!;
    expect(within(modal).getByRole('link', { name: /^GenomeArk/ })).toHaveAttribute(
      'href',
      'https://www.genomeark.org/genomeark-all/Acipenser_ruthenus.html',
    );
    await userEvent.click(screen.getByRole('button', { name: 'Close' }));
    expect(
      screen.queryByText('No quality metrics have been published for this assembly.'),
    ).not.toBeInTheDocument();
  });

  it('never claims to be official or endorsed', async () => {
    seed();
    const { container } = renderAt('/collections/vgp');
    await loaded();
    expect(container.textContent).not.toMatch(/official/i);
    expect(container.textContent).not.toMatch(/endors/i);
  });

  it('says so for an unknown collection, with a way back', () => {
    renderAt('/collections/unknown');
    expect(screen.getByText(/No collection named “unknown”/)).toBeInTheDocument();
    expect(screen.getByRole('link', { name: 'See all collections' })).toHaveAttribute(
      'href',
      '/collections',
    );
  });
});
