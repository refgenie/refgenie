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
  { digest: OTHER, label: 'rCRSd · GenBank' },
  { digest: THIRD, label: 'hg38_primary · UCSC' },
];

describe('SequencesPanel', () => {
  it('names the alphabet and shows the table', async () => {
    renderWithProviders(<SequencesPanel digest={DIGEST} />);
    const figure = (await screen.findByText('Alphabets')).closest('figure');
    expect(figure).toHaveTextContent('dna3bit');
    expect(figure).toHaveTextContent('dna2bit');
    expect(screen.getByRole('table', { name: 'Sequences' })).toBeInTheDocument();
    expect(screen.getByText(/SHA-512\/24u, the GA4GH refget ID/)).toBeInTheDocument();
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

function renderExtract(route = '/?tab=extract', browser = true) {
  return renderWithProviders(
    <ExtractPanel digest={DIGEST} storeUrl={STORE} browserStoreUrl={browser ? STORE : null} />,
    { route },
  );
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
        <ExtractPanel digest={DIGEST} storeUrl={STORE} browserStoreUrl={STORE} />
      </StrictMode>,
      { route: '/?region=extra:1-4' },
    );
    expect(await screen.findByText(/^>extra:1-4/)).toBeInTheDocument();
  });

  it('says the browser reads up to 100,000 bases', async () => {
    renderExtract();
    expect(await screen.findByText(/Up to 100,000 bases/)).toBeInTheDocument();
  });

  it('shows code for the region typed, without the browser box when it is off', async () => {
    renderExtract('/?tab=extract&region=extra:1-4', false);
    expect(await screen.findByRole('heading', { name: 'In your own code' })).toBeInTheDocument();
    expect(screen.queryByLabelText('Region')).toBeNull();
    expect(screen.getByText(/--start 0 --end 4 --remote https:\/\/store\.test\/s\/$/)).toBeInTheDocument();
  });

  it('keeps the code in step with the region box', async () => {
    renderExtract('/?tab=extract&region=extra:1-4');
    await screen.findByText(/^>extra:1-4/);
    const input = screen.getByLabelText('Region');
    await userEvent.clear(input);
    await userEvent.type(input, 'rCRSd:11-20');
    expect(screen.getByText(/--start 10 --end 20 /)).toBeInTheDocument();
    // Half-typed text keeps the last region that parsed.
    await userEvent.type(input, ':');
    expect(screen.getByText(/--start 10 --end 20 /)).toBeInTheDocument();
    await userEvent.click(screen.getByRole('tab', { name: 'Python' }));
    expect(screen.getByText(/get_substring\(".*", 10, 20\)/)).toBeInTheDocument();
  });
});

async function bedMode(route = '/?tab=extract', browser = true) {
  renderExtract(route, browser);
  await userEvent.click(await screen.findByRole('radio', { name: 'BED file' }));
}

function pasteBed(text: string) {
  fireEvent.change(screen.getByLabelText('BED lines'), { target: { value: text } });
}

describe('ExtractPanel BED mode', () => {
  it('switches between one region and a BED file', async () => {
    renderExtract();
    const one = await screen.findByRole('radio', { name: 'One region' });
    expect(one).toHaveAttribute('aria-checked', 'true');
    await userEvent.click(screen.getByRole('radio', { name: 'BED file' }));
    expect(screen.getByRole('radio', { name: 'BED file' })).toHaveAttribute(
      'aria-checked',
      'true',
    );
    expect(screen.getByLabelText('BED lines')).toBeVisible();
    expect(screen.getByLabelText('Region')).not.toBeVisible();
    expect(screen.getByText(/Up to 5,000 regions and 5,000,000 bases/)).toBeInTheDocument();
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
    expect(within(problems).getByText(/Line 4/).closest('li')).toHaveTextContent(
      'Line 4 Extra 0 4 No sequence named "Extra". Did you mean "extra"?',
    );
    expect(problems).toHaveTextContent('Line 5');
    expect(problems).toHaveTextContent('past the end of extra (80 bases)');

    await userEvent.click(screen.getByRole('button', { name: 'Extract 2 regions' }));
    expect(await screen.findByText(/^>peak::NC_012920.1:1-9/)).toHaveTextContent(
      '>peak::NC_012920.1:1-9 CGTACGTA >extra:0-4 ACGT',
    );
    expect(screen.getByText(/2 records · 12 bases/)).toBeInTheDocument();
    expect(screen.getByRole('button', { name: 'Copy the FASTA' })).toBeInTheDocument();
    expect(screen.getByRole('button', { name: /Download .fa/ })).toBeInTheDocument();
  });

  it('reads a dropped BED file', async () => {
    await bedMode();
    const file = new File(['extra\t0\t8\tone\nrCRSd\t4\t8\n'], 'peaks.bed', {
      type: 'text/plain',
    });
    fireEvent.drop(screen.getByTestId('bed-drop'), { dataTransfer: { files: [file] } });
    expect(await screen.findByText('peaks.bed')).toBeInTheDocument();
    expect(screen.getByText('2 regions')).toBeInTheDocument();
    await userEvent.click(screen.getByRole('button', { name: 'Extract 2 regions' }));
    expect(await screen.findByText(/^>one::extra:0-8/)).toBeInTheDocument();
  });

  it('reads a chosen file, and can remove it', async () => {
    await bedMode();
    const file = new File(['extra\t0\t8\n'], 'mine.bed');
    await userEvent.upload(screen.getByLabelText('Choose a BED file'), file);
    expect(await screen.findByText('mine.bed')).toBeInTheDocument();
    expect(screen.getByText('1 region')).toBeInTheDocument();
    await userEvent.click(screen.getByRole('button', { name: 'Remove' }));
    expect(screen.getByLabelText('BED lines')).toBeInTheDocument();
  });

  it('fills in an example from this genome', async () => {
    await bedMode();
    await userEvent.click(screen.getByRole('button', { name: 'Try an example' }));
    expect(screen.getByLabelText('BED lines')).toHaveValue(
      'rCRSd\t16318\t16819\tregion_1\nextra\t40\t80\tregion_2\nextra\t0\t60\n',
    );
    expect(screen.getByRole('button', { name: 'Extract 3 regions' })).toBeEnabled();
  });

  it('will not run past the region limit', async () => {
    await bedMode();
    pasteBed(Array.from({ length: 5001 }, () => 'extra\t0\t1').join('\n'));
    expect(screen.getByRole('alert')).toHaveTextContent(
      'This file has 5,001 regions (the limit is 5,000).',
    );
    expect(screen.getByRole('button', { name: 'Extract 5,001 regions' })).toBeDisabled();
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

  it('shows BED code for each tool, without the browser box when it is off', async () => {
    await bedMode('/?tab=extract', false);
    expect(screen.queryByLabelText('BED lines')).toBeNull();
    expect(
      screen.getByText(`refget store export ${DIGEST} --bed regions.bed --remote ${STORE} -o regions.fa`),
    ).toBeInTheDocument();
    await userEvent.click(screen.getByRole('tab', { name: 'Python' }));
    expect(screen.getByText(/export_fasta_from_regions/)).toBeInTheDocument();
  });
});

describe('ComparePanel', () => {
  it('does not show the collection digests', async () => {
    renderWithProviders(<ComparePanel digest={DIGEST} candidates={CANDIDATES} />);
    await screen.findByLabelText(/Other genome/);
    expect(screen.queryByRole('heading', { name: 'Collection digests' })).toBeNull();
  });

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

  it('shows the server summary in place of its own sentences, and lists renames', async () => {
    const summary =
      'All 1 sequences of hg38 are in rCRSd, which has 1 more; 1 names differ: chrM->MT.';
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
