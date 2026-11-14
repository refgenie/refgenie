import { describe, expect, it, onTestFinished, vi } from 'vitest';
import { screen, waitFor, within } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { HttpResponse, http } from 'msw';
import { Route, Routes } from 'react-router-dom';
import { GenomePage } from './GenomePage';
import { server } from '../test/server';
import { API } from '../test/handlers';
import { renderWithProviders, serverConfig } from '../test/renderWithProviders';
import { emptyPage, genomeDetailFixture, genomesFixture } from '../test/fixtures';

const flags = vi.hoisted(() => ({ fasta: true }));
vi.mock('../services/config', async (importOriginal) => ({
  ...(await importOriginal<typeof import('../services/config')>()),
  fastaDownloadEnabled: () => flags.fasta,
}));

vi.mock('../services/fastaExport', async (importOriginal) => ({
  ...(await importOriginal<typeof import('../services/fastaExport')>()),
  loadCollection: async () => ({ seqTemplate: 'sequences/%s2/%s.seq', seqs: [] }),
}));

const DIGEST = '5cb5b9d2d1e2b5e5e2f0a3c4d5e6f708';

function renderGenome(route = `/genomes/${DIGEST}`) {
  return renderWithProviders(
    <Routes>
      <Route path="/genomes/:digest" element={<GenomePage />} />
    </Routes>,
    { route },
  );
}

describe('GenomePage', () => {
  it('renders the genome detail parity fields on the Details tab', async () => {
    renderGenome(`/genomes/${DIGEST}?tab=details`);
    expect(
      await screen.findByRole('heading', { name: 'hg38', level: 1 }),
    ).toBeInTheDocument();
    expect(screen.getByText('chromosome')).toBeInTheDocument();
    // Taxon and accession are rendered as external links. The name carries the
    // sr-only "(opens in a new tab)" suffix that ExternalLink appends, because
    // the arrow icon itself is aria-hidden.
    expect(screen.getByRole('link', { name: /^9606\b/ })).toHaveAttribute(
      'href',
      'https://identifiers.org/taxonomy:9606',
    );
    const accessions = screen.getAllByRole('link', { name: /^GCA_000001405\.15\b/ });
    // Once in the head, once in Details.
    expect(accessions).toHaveLength(2);
    expect(accessions[1]).toHaveAttribute(
      'href',
      'https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_000001405.15/',
    );
  });

  it("shows the species list's other names", async () => {
    renderGenome(`/genomes/${DIGEST}?tab=details`);
    await screen.findByRole('heading', { name: 'hg38', level: 1 });
    const term = screen.getByText('Other species names');
    expect(term.nextElementSibling).toHaveTextContent('modern human');
  });

  it('shows the metadata panel on the Details tab', async () => {
    renderGenome(`/genomes/${DIGEST}?tab=details`);
    expect(
      await screen.findByRole('heading', { name: 'Metadata' }),
    ).toBeInTheDocument();
  });

  it('opens the raw metadata JSON in a modal', async () => {
    renderGenome(`/genomes/${DIGEST}?tab=details`);
    await userEvent.click(await screen.findByRole('button', { name: 'Raw metadata JSON' }));
    const dialog = screen.getByRole('dialog', { name: 'Raw metadata JSON' });
    expect(dialog).toHaveTextContent('"license": "CC0-1.0"');
  });

  it('groups every digest in one refget section on the Details tab', async () => {
    renderGenome(`/genomes/${DIGEST}?tab=details`);
    const heading = await screen.findByRole('heading', { name: /refget sequence collection/ });
    const section = heading.closest('section') as HTMLElement;
    expect(within(section).getByText('Top-level collection digest')).toBeInTheDocument();
    expect(within(section).getByText('names')).toBeInTheDocument();
    expect(within(section).getByRole('button', { name: 'Copy names digest' })).toBeInTheDocument();
    expect(
      within(section).getAllByRole('button', { name: 'Genomes that share this' }).length,
    ).toBeGreaterThan(0);
  });

  it('lists file-mode assets that have no archive record', async () => {
    renderGenome(`/genomes/${DIGEST}?tab=assets`);
    // fasta_index is serving_modes: ["file"] and has no tarball, but the table
    // is driven from /assets so it still gets a row.
    expect(await screen.findByText('FASTA index')).toBeInTheDocument();
  });

  it('offers the FASTA download on the Download tab, not in the head', async () => {
    renderGenome();
    await screen.findByRole('tab', { name: 'Download', selected: false });
    expect(screen.queryByRole('button', { name: /^Download/ })).toBeNull();
    await userEvent.click(screen.getByRole('tab', { name: 'Download' }));
    expect(
      await screen.findByRole('button', { name: /^Download \.fa\.gz/ }),
    ).toBeInTheDocument();
    expect(screen.getByRole('radio', { name: 'All sequences' })).toBeChecked();
    expect(screen.getByRole('tab', { name: 'refget CLI' })).toBeInTheDocument();
  });

  it('offers no FASTA download without a refget store', async () => {
    server.use(
      http.get(`${API}/genomes/:digest`, () =>
        HttpResponse.json({ ...genomeDetailFixture, refget_store_url: null }),
      ),
    );
    renderGenome(`/genomes/${DIGEST}?tab=download`);
    await screen.findByRole('heading', { name: 'hg38', level: 1 });
    expect(screen.queryByRole('button', { name: /^Download/ })).toBeNull();
    for (const name of ['Download', 'Extract']) {
      expect(screen.getByRole('tab', { name })).toBeDisabled();
      expect(screen.getByRole('tab', { name })).toHaveAttribute(
        'title',
        'This genome has no refget store',
      );
    }
  });

  it('opens the Extract tab with the region in the URL', async () => {
    renderGenome(`/genomes/${DIGEST}?tab=extract&region=extra:1-4`);
    expect(await screen.findByRole('tab', { selected: true })).toHaveTextContent(
      'Extract',
    );
    expect(await screen.findByLabelText('Region')).toHaveValue('extra:1-4');
    expect(screen.getByRole('heading', { name: 'In your own code' })).toBeInTheDocument();
  });

  it('shows the whole description, with no More toggle', async () => {
    renderGenome();
    await screen.findByRole('heading', { name: 'hg38', level: 1 });
    expect(screen.queryByRole('button', { name: /^(More|Less)$/ })).toBeNull();
  });

  it('shows key facts in the head', async () => {
    renderGenome();
    await screen.findByRole('heading', { name: 'hg38', level: 1 });
    expect(await screen.findByText('refget digest')).toBeInTheDocument();
    expect(screen.queryByText(/built assets?$/)).toBeNull();
    expect(screen.getByRole('link', { name: /What is a refget digest/ })).toHaveAttribute(
      'href',
      'https://ga4gh.github.io/refget/seqcols/',
    );
  });

  it('shows the size of the fasta pull beside its command', async () => {
    renderWithProviders(
      <Routes>
        <Route path="/genomes/:digest" element={<GenomePage />} />
      </Routes>,
      { route: `/genomes/${DIGEST}`, config: serverConfig },
    );
    expect(await screen.findByText('refgenie pull hg38/fasta')).toBeInTheDocument();
    // 3,273,481,150 bytes in the assets fixture.
    expect(screen.getByTitle('Size of the asset once pulled')).toHaveTextContent('3.0 GB');
  });

  it('keeps the alias list off the head; Details has it, each with a copy button', async () => {
    const names = ['hg38', 'GRCh38-a', 'GRCh38-b', 'hs38', 'hg38-x', 'hg38-y'];
    server.use(
      http.get(`${API}/aliases`, () =>
        HttpResponse.json({
          results: names.map((name) => ({ name, genome_digest: DIGEST })),
          pagination: { page: 0, page_size: 100, total: names.length },
        }),
      ),
    );
    renderGenome();
    await screen.findByRole('heading', { name: 'hg38', level: 1 });
    expect(screen.queryByText(/Also known as/)).toBeNull();
    expect(screen.queryByRole('link', { name: /more$/ })).toBeNull();
    await userEvent.click(screen.getByRole('tab', { name: 'Details' }));
    const table = await screen.findByRole('table', { name: 'Aliases' });
    expect(within(table).getAllByRole('row')).toHaveLength(names.length + 2);
    // The button sits in the alias's own cell.
    const cell = within(table).getByText('hs38').closest('td') as HTMLElement;
    expect(within(cell).getByRole('button', { name: 'Copy hs38' })).toBeInTheDocument();
  });

  it('lists genome and store aliases on the Details tab, with their source', async () => {
    renderGenome(`/genomes/${DIGEST}?tab=details`);
    const table = await screen.findByRole('table', { name: 'Aliases' });
    expect(within(table).getByText('hg38').closest('tr')).toHaveTextContent('refgenie');
    expect(
      (await within(table).findByText('rCRSd_extra')).closest('tr'),
    ).toHaveTextContent('refget store, refgenie');
  });

  it('renders 404 copy naming the digest', async () => {
    server.use(
      http.get(`${API}/genomes/:digest`, () =>
        HttpResponse.json({ detail: 'Genome not found' }, { status: 404 }),
      ),
    );
    renderGenome();
    const alert = await screen.findByRole('alert');
    expect(alert).toHaveTextContent('Not found');
    expect(alert).toHaveTextContent(DIGEST);
  });

  // The head is a sibling of the loading/error branches, not something they
  // return past. A missing genome must still say where you are and offer a way
  // back, and must never title the tab with the "Genome" placeholder.
  it('keeps the head, and the way out, when the genome is missing', async () => {
    server.use(
      http.get(`${API}/genomes/:digest`, () =>
        HttpResponse.json({ detail: 'Genome not found' }, { status: 404 }),
      ),
    );
    renderGenome();
    await screen.findByRole('alert');
    expect(screen.getByRole('heading', { level: 1 })).toHaveTextContent('Genome');
    expect(screen.getByRole('link', { name: 'Genomes' })).toHaveAttribute(
      'href',
      '/genomes',
    );
    expect(document.title).toBe('refgenie local dashboard');
  });

  // A slow remote-store load used to time out at the gateway and leave the page
  // on "Loading genome…". The error must show, with a way to try again.
  it('shows an error with a retry action on a 504', async () => {
    let calls = 0;
    server.use(
      http.get(`${API}/genomes/:digest`, () => {
        calls += 1;
        return calls === 1
          ? HttpResponse.text('Gateway Time-out', { status: 504 })
          : HttpResponse.json({ detail: 'Genome not found' }, { status: 404 });
      }),
    );
    renderGenome();
    const alert = await screen.findByRole('alert');
    expect(alert).toHaveTextContent('504');
    expect(screen.queryByText(/Loading genome/)).not.toBeInTheDocument();
    await userEvent.click(screen.getByRole('button', { name: 'Retry' }));
    expect(await screen.findByText(/Not found/)).toBeInTheDocument();
    expect(calls).toBe(2);
  });

  it('opens on the Assets tab', async () => {
    renderGenome();
    const tabs = await screen.findAllByRole('tab');
    expect(tabs.map((tab) => tab.textContent?.replace(/\d|,/g, ''))).toEqual([
      'Sequences',
      'Extract',
      'Download',
      'Assets',
      'Compare',
      'Details',
    ]);
    expect(await screen.findByRole('tab', { selected: true })).toHaveTextContent(
      /^Assets/,
    );
  });

  it('opens the tab named in the URL', async () => {
    renderGenome(`/genomes/${DIGEST}?tab=sequences`);
    expect(await screen.findByRole('tab', { selected: true })).toHaveTextContent(
      /^Sequences/,
    );
    expect(await screen.findByRole('table', { name: 'Sequences' })).toBeInTheDocument();
  });

  it('moves between tabs with the arrow keys', async () => {
    renderGenome(`/genomes/${DIGEST}?tab=details`);
    const details = await screen.findByRole('tab', { name: 'Details' });
    details.focus();
    await userEvent.keyboard('{ArrowRight}');
    expect(screen.getByRole('tab', { selected: true })).toHaveTextContent(/^Sequences/);
    expect(document.activeElement).toBe(screen.getByRole('tab', { selected: true }));
    await userEvent.keyboard('{ArrowLeft}{ArrowLeft}');
    expect(screen.getByRole('tab', { selected: true })).toHaveTextContent(/^Compare/);
  });

  it('disables Sequences and Compare when the genome has no sequence data', async () => {
    server.use(
      http.get(
        `${API}/genomes/:digest/sequences`,
        () => new HttpResponse(null, { status: 404 }),
      ),
    );
    renderGenome(`/genomes/${DIGEST}?tab=sequences`);
    await waitFor(() =>
      expect(screen.getByRole('tab', { name: /^Sequences/ })).toBeDisabled(),
    );
    expect(screen.getByRole('tab', { name: 'Compare' })).toBeDisabled();
    // The URL asked for a tab that has nothing, so the default opens instead.
    expect(screen.getByRole('tab', { selected: true })).toHaveTextContent(/^Assets/);
    expect(screen.queryByText(/not available/)).not.toBeInTheDocument();
  });

  it('asks for genomes of the same taxon, and counts them on the Compare tab', async () => {
    let query: URLSearchParams | undefined;
    server.use(
      http.get(`${API}/genomes`, ({ request }) => {
        query = new URL(request.url).searchParams;
        return HttpResponse.json(genomesFixture);
      }),
    );
    renderGenome(`/genomes/${DIGEST}?tab=compare`);
    const compare = await screen.findByRole('tab', { name: 'Compare 1' });
    expect(compare).toBeEnabled();
    expect(query?.get('taxon_id')).toBe('9606');
  });

  it('disables Compare when no other genome shares the species', async () => {
    flags.fasta = false;
    onTestFinished(() => {
      flags.fasta = true;
    });
    server.use(
      http.get(`${API}/genomes`, () =>
        HttpResponse.json({
          ...genomesFixture,
          results: genomesFixture.results.filter((g) => g.digest === DIGEST),
        }),
      ),
    );
    renderGenome(`/genomes/${DIGEST}?tab=compare`);
    await waitFor(() =>
      expect(screen.getByRole('tab', { name: 'Compare' })).toBeDisabled(),
    );
    expect(screen.getByRole('tab', { name: 'Compare' })).toHaveAttribute(
      'title',
      'No other genomes of this species',
    );
  });

  it('keeps Compare open for your own FASTA when no other genome shares the species', async () => {
    server.use(
      http.get(`${API}/genomes`, () =>
        HttpResponse.json({
          ...genomesFixture,
          results: genomesFixture.results.filter((g) => g.digest === DIGEST),
        }),
      ),
    );
    renderGenome(`/genomes/${DIGEST}?tab=compare`);
    expect(
      await screen.findByRole('heading', { name: 'Compare your own FASTA' }),
    ).toBeInTheDocument();
    expect(screen.getByRole('tab', { name: 'Compare' })).toBeEnabled();
    expect(
      await screen.findByText('No other genomes of this species to compare with.'),
    ).toBeInTheDocument();
  });

  it('on a server with no assets, shows how to request them', async () => {
    server.use(http.get(`${API}/assets`, () => HttpResponse.json(emptyPage)));
    renderWithProviders(
      <Routes>
        <Route path="/genomes/:digest" element={<GenomePage />} />
      </Routes>,
      { route: `/genomes/${DIGEST}`, config: serverConfig },
    );
    expect(await screen.findByText('No assets yet')).toBeInTheDocument();
    const link = screen.getByRole('link', { name: /Request assets/ });
    const href = link.getAttribute('href') ?? '';
    expect(href).toContain('genome=');
    expect(href).toContain('genome_digest=');
    expect(screen.queryByText('Build the first one')).not.toBeInTheDocument();
  });

  it('on a local dash with no assets, offers to build and no request link', async () => {
    server.use(http.get(`${API}/assets`, () => HttpResponse.json(emptyPage)));
    renderGenome();
    expect(await screen.findByText('Build the first one')).toBeInTheDocument();
    expect(screen.queryByRole('link', { name: /Request assets/ })).not.toBeInTheDocument();
  });
});
