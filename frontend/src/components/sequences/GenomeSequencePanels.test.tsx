import { File as NodeFile } from 'node:buffer';
import { StrictMode } from 'react';
import { describe, expect, it, vi } from 'vitest';
import { fireEvent, screen, within } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { HttpResponse, http } from 'msw';
import { ComparePanel, ExtractPanel, SequencesPanel } from './GenomeSequencePanels';
import { server } from '../../test/server';
import { API } from '../../test/handlers';
import { renderWithProviders } from '../../test/renderWithProviders';
import { genomeComparisonFixture } from '../../test/fixtures';

const extractRegion = vi.hoisted(() =>
  vi.fn(async ({ start, end }: { start: number; end: number }) =>
    'ACGT'.repeat(50).slice(start, end),
  ),
);
vi.mock('../../services/regionExtract', async (importOriginal) => ({
  ...(await importOriginal<typeof import('../../services/regionExtract')>()),
  extractRegion,
}));

const STORE = 'https://store.test/s/';
const DIGEST = '5cb5b9d2d1e2b5e5e2f0a3c4d5e6f708';
const OTHER = 'a1b2c3d4e5f60718293a4b5c6d7e8f90';
const THIRD = 'ffffffffffffffffffffffffffffffff';
const CANDIDATES = [
  { digest: OTHER, label: 'rCRSd · GenBank', aliases: ['rCRSd', 'mito'] },
  { digest: THIRD, label: 'hg38_primary · UCSC', aliases: ['hg38_primary'] },
];

describe('SequencesPanel', () => {
  it('names the alphabet and shows the table', async () => {
    renderWithProviders(<SequencesPanel digest={DIGEST} />);
    const figure = (await screen.findByText('Alphabets')).closest('figure');
    expect(figure).toHaveTextContent('dna3bit');
    expect(figure).toHaveTextContent('dna2bit');
    expect(screen.getByRole('table', { name: 'Sequences' })).toBeInTheDocument();
  });
});

describe('SequencesPanel extract links', () => {
  it('has none unless the genome can be extracted', async () => {
    renderWithProviders(<SequencesPanel digest={DIGEST} />);
    await screen.findByRole('table', { name: 'Sequences' });
    expect(screen.queryByRole('link', { name: /Extract a region/ })).toBeNull();
  });

  it('links each row to the Extract tab with its first bases', async () => {
    renderWithProviders(<SequencesPanel digest={DIGEST} canExtract />);
    expect(
      await screen.findByRole('link', { name: 'Extract a region of extra' }),
    ).toHaveAttribute('href', '/?tab=extract&region=extra%3A1-80');
  });
});

function renderExtract(route = '/?tab=extract') {
  return renderWithProviders(<ExtractPanel digest={DIGEST} storeUrl={STORE} />, { route });
}

describe('ExtractPanel', () => {
  it('extracts a region as FASTA, matching names by alias', async () => {
    renderExtract();
    const input = await screen.findByLabelText('Region');
    await userEvent.type(input, 'NC_012920.1:2-9');
    await userEvent.click(screen.getByRole('button', { name: 'Extract' }));
    expect(await screen.findByText(/^>rCRSd:2-9/)).toHaveTextContent(
      '>rCRSd:2-9 CGTACGTA',
    );
    expect(extractRegion).toHaveBeenLastCalledWith(
      expect.objectContaining({ baseUrl: STORE, start: 1, end: 9 }),
    );
    expect(screen.getByText(/8 bases/)).toBeInTheDocument();
    expect(screen.getByRole('button', { name: 'Copy the FASTA' })).toBeInTheDocument();
  });

  it('shows an error with a suggestion for an unknown name', async () => {
    renderExtract();
    await userEvent.type(await screen.findByLabelText('Region'), 'rcrsd:1-10{Enter}');
    expect(await screen.findByRole('alert')).toHaveTextContent('Did you mean "rCRSd"?');
  });

  it('runs a region from the URL', async () => {
    renderExtract('/?tab=extract&region=extra:1-4');
    expect(await screen.findByText(/^>extra:1-4/)).toHaveTextContent('>extra:1-4 ACGT');
    expect(screen.getByLabelText('Region')).toHaveValue('extra:1-4');
  });

  it('runs a region from the URL under StrictMode', async () => {
    renderWithProviders(
      <StrictMode>
        <ExtractPanel digest={DIGEST} storeUrl={STORE} />
      </StrictMode>,
      { route: '/?region=extra:1-4' },
    );
    expect(await screen.findByText(/^>extra:1-4/)).toBeInTheDocument();
  });

  it('keeps the code in step with the region box', async () => {
    renderExtract('/?tab=extract&region=extra:1-4');
    await screen.findByText(/^>extra:1-4/);
    const input = screen.getByLabelText('Region');
    await userEvent.clear(input);
    await userEvent.type(input, 'rCRSd:11-20');
    expect(screen.getByText(/--start 10 --end 20 /)).toBeInTheDocument();
  });
});

async function bedMode(route = '/?tab=extract') {
  renderExtract(route);
  await userEvent.click(await screen.findByRole('radio', { name: 'BED file' }));
}

// jsdom's File has no arrayBuffer(); Node's does, like a browser's.
const bedFile = (text: string, name: string) => new NodeFile([text], name) as unknown as File;

function pasteBed(text: string) {
  fireEvent.change(screen.getByLabelText('Or paste BED lines'), {
    target: { value: text },
  });
}

describe('ExtractPanel BED mode', () => {
  it('switches between one region and a BED file', async () => {
    renderExtract();
    expect(await screen.findByRole('radio', { name: 'One region' })).toBeChecked();
    await userEvent.click(screen.getByRole('radio', { name: 'BED file' }));
    expect(screen.getByRole('radio', { name: 'BED file' })).toBeChecked();
    expect(screen.getByLabelText('Or paste BED lines')).toBeVisible();
    expect(screen.getByLabelText('Region')).not.toBeVisible();
  });

  it('summarizes a pasted BED, lists problem lines, and runs the good ones', async () => {
    await bedMode();
    pasteBed(
      ['track name=x', 'NC_012920.1\t1\t9\tpeak', 'extra 0 4', 'Extra\t0\t4', 'extra\t70\t90'].join(
        '\n',
      ),
    );
    expect(screen.getByText('2 regions')).toBeInTheDocument();
    expect(screen.getByText(/12 bases/)).toBeInTheDocument();
    const problems = screen.getByRole('list', { name: 'Problem lines' });
    expect(within(problems).getAllByRole('listitem')).toHaveLength(2);

    await userEvent.click(screen.getByRole('button', { name: 'Extract 2 regions' }));
    expect(await screen.findByText(/^>peak::NC_012920.1:1-9/)).toHaveTextContent(
      '>peak::NC_012920.1:1-9 CGTACGTA >extra:0-4 ACGT',
    );
    expect(screen.getByText(/2 records · 12 bases/)).toBeInTheDocument();
    expect(screen.getByRole('button', { name: 'Copy the FASTA' })).toBeInTheDocument();
    expect(screen.getByRole('button', { name: /Download .fa/ })).toBeInTheDocument();
  });

  it('reads a chosen file, and can remove it', async () => {
    await bedMode();
    const file = bedFile('extra\t0\t8\tone\nrCRSd\t4\t8\n', 'peaks.bed');
    await userEvent.upload(screen.getByLabelText(/^BED file \(/), file);
    expect(await screen.findByText('2 regions')).toBeInTheDocument();
    expect(screen.queryByLabelText('Or paste BED lines')).toBeNull();
    await userEvent.click(screen.getByRole('button', { name: 'Extract 2 regions' }));
    expect(await screen.findByText(/^>one::extra:0-8/)).toBeInTheDocument();
    await userEvent.click(screen.getByRole('button', { name: 'Remove file' }));
    expect(screen.getByLabelText('Or paste BED lines')).toBeInTheDocument();
  });

  it('shows progress and can be cancelled', async () => {
    extractRegion.mockImplementationOnce(
      (args: { start: number; end: number; signal?: AbortSignal }) =>
        new Promise<string>((_, reject) => {
          args.signal?.addEventListener('abort', () => reject(new Error('aborted')));
        }),
    );
    await bedMode();
    pasteBed('extra\t0\t4\nextra\t4\t8\n');
    await userEvent.click(screen.getByRole('button', { name: 'Extract 2 regions' }));
    expect(await screen.findByRole('progressbar', { name: 'Regions read' })).toBeInTheDocument();
    await userEvent.click(screen.getByRole('button', { name: 'Cancel' }));
    expect(screen.getByText('Cancelled.')).toBeInTheDocument();
    expect(screen.queryByText(/^>extra/)).toBeNull();
  });
});

describe('ComparePanel', () => {
  it('compares with another genome', async () => {
    let requested: string | undefined;
    server.use(
      http.get(`${API}/genomes/:digest/compare/:other`, ({ params }) => {
        requested = `${params.digest}/${params.other}`;
        return HttpResponse.json(genomeComparisonFixture);
      }),
    );
    renderWithProviders(<ComparePanel digest={DIGEST} candidates={CANDIDATES} />);
    const box = await screen.findByRole('combobox', { name: /Other genome/ });
    await userEvent.type(box, 'rCRSd');
    expect(requested).toBeUndefined();
    await userEvent.type(box, ' · GenBank');
    expect(
      await screen.findByText(
        "This genome's coordinate system is a subset of the other genome's.",
      ),
    ).toBeInTheDocument();
    expect(requested).toBe(`${DIGEST}/${OTHER}`);
    expect(screen.getByRole('link', { name: 'rCRSd · GenBank' })).toHaveAttribute(
      'href',
      `/genomes/${OTHER}`,
    );
  });

  it('compares with the genome an alias or partial name picks on Enter', async () => {
    let requested: string | undefined;
    server.use(
      http.get(`${API}/genomes/:digest/compare/:other`, ({ params }) => {
        requested = String(params.other);
        return HttpResponse.json(genomeComparisonFixture);
      }),
    );
    renderWithProviders(<ComparePanel digest={DIGEST} candidates={CANDIDATES} />);
    const box = await screen.findByRole('combobox', { name: /Other genome/ });
    await userEvent.type(box, 'MITO{Enter}');
    expect(
      await screen.findByRole('link', { name: 'rCRSd · GenBank' }),
    ).toBeInTheDocument();
    expect(requested).toBe(OTHER);
    await userEvent.clear(box);
    await userEvent.type(box, 'primary{Enter}');
    expect(box).toHaveValue('hg38_primary · UCSC');
    expect(
      await screen.findByRole('link', { name: 'hg38_primary · UCSC' }),
    ).toBeInTheDocument();
  });

  it('says so when nothing is picked', async () => {
    renderWithProviders(<ComparePanel digest={DIGEST} candidates={CANDIDATES} />);
    await userEvent.click(await screen.findByRole('button', { name: 'Compare' }));
    expect(screen.getByRole('alert')).toHaveTextContent(
      'No genome selected. Choose one from the list.',
    );
    const box = screen.getByRole('combobox', { name: /Other genome/ });
    await userEvent.type(box, 'mm10{Enter}');
    expect(screen.getByRole('alert')).toHaveTextContent(
      'No genome of this species matches "mm10".',
    );
  });

  it('shows the server summary in place of its own sentences, and lists renames', async () => {
    const summary =
      'The one sequence of hg38 is in rCRSd, which has 1 more; 1 name differs: chrM->MT.';
    server.use(
      http.get(`${API}/genomes/:digest/compare/:other`, () =>
        HttpResponse.json({
          ...genomeComparisonFixture,
          summary,
          name_map: { chrM: 'MT' },
        }),
      ),
    );
    renderWithProviders(<ComparePanel digest={DIGEST} candidates={CANDIDATES} />);
    await userEvent.type(
      await screen.findByRole('combobox', { name: /Other genome/ }),
      'rCRSd · GenBank',
    );
    expect(await screen.findByText(summary)).toBeInTheDocument();
    expect(screen.queryByText(/coordinate system is a subset/)).toBeNull();
    const renames = screen.getByRole('table', { name: '1 sequence has another name' });
    expect(within(renames).getByText('MT')).toBeInTheDocument();
  });
});
