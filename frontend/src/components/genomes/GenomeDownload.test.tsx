import { beforeEach, describe, expect, it, vi } from 'vitest';
import { screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { HttpResponse, http } from 'msw';
import { useState } from 'react';
import { GenomeDownloadPanel } from './GenomeDownload';
import { useFastaDownload } from '../../hooks/useFastaDownload';
import { genomeDetailFixture, genomeSequencesFixture } from '../../test/fixtures';
import { renderWithProviders as render } from '../../test/renderWithProviders';
import { server } from '../../test/server';
import { API } from '../../test/handlers';
import type { GenomeDetailResponse } from '../../types/api';

vi.mock('../../services/saveStream', () => ({
  // Drain the chunks, as a real save would.
  saveStream: vi.fn(
    async ({ chunks }: { chunks: () => AsyncGenerator<Uint8Array> }) => {
      for await (const chunk of chunks()) void chunk;
    },
  ),
}));

vi.mock('../../services/fastaExport', async (importOriginal) => ({
  ...(await importOriginal<typeof import('../../services/fastaExport')>()),
  loadGtars: vi.fn(async () => ({})),
  fastaChunks: vi.fn(async function* () {
    yield new TextEncoder().encode('>chr1\nACGT\n');
  }),
  loadSeqTemplate: vi.fn(async () => 'sequences/%s2/%s.seq'),
}));

const { saveStream } = await import('../../services/saveStream');
const { fastaChunks } = await import('../../services/fastaExport');

const row = (name: string, alphabet = 'dna2bit') => ({
  name,
  length: 80,
  alphabet,
  description: null,
  sha512t24u: name,
  md5: null,
  aliases: [],
});

/** The genome's `/sequences` rows. */
function serveSequences(rows: ReturnType<typeof row>[]) {
  server.use(
    http.get(`${API}/genomes/:digest/sequences`, () =>
      HttpResponse.json({ ...genomeSequencesFixture, sequences: rows }),
    ),
  );
}

const PICK_SEQS = ['chr1', 'chr2', 'chr1_KI270706v1_random', 'chrEBV'].map((name) =>
  row(name),
);

const genome: GenomeDetailResponse = genomeDetailFixture;

function Harness({ storeUrl }: { storeUrl: string | null }) {
  const download = useFastaDownload({
    storeUrl,
    digest: genome.digest,
    fileStem: 'hg38',
  });
  const [open, setOpen] = useState(true);
  return (
    <>
      <button type="button" onClick={() => setOpen(!open)}>
        Toggle
      </button>
      {open && (
        <GenomeDownloadPanel
          digest={genome.digest}
          storeUrl={storeUrl ?? ''}
          alias="hg38"
          download={download}
          hasFastaAsset
        />
      )}
    </>
  );
}

const renderPanel = () => render(<Harness storeUrl={genome.refget_store_url} />);

describe('Download tab', () => {
  beforeEach(() => {
    vi.mocked(saveStream).mockClear();
    serveSequences([row('chr1')]);
  });

  it('saves to <fileStem>.fa.gz', async () => {
    renderPanel();
    expect(screen.getByRole('heading', { name: 'In your browser' })).toBeInTheDocument();
    // Enabled once the collection has loaded and the size is known.
    expect(await screen.findByText(/\(87 B uncompressed\)/)).toBeInTheDocument();
    await userEvent.click(screen.getByRole('button', { name: /^Download \.fa\.gz/ }));
    expect(saveStream).toHaveBeenCalledWith(
      expect.objectContaining({ suggestedName: 'hg38.fa.gz', estimatedBytes: 87 }),
    );
    expect(await screen.findByText('Saved')).toBeInTheDocument();
  });

  it('exports only the chosen sequences, and the CLI command names them', async () => {
    serveSequences(PICK_SEQS);
    renderPanel();
    await screen.findByText(/4 sequences/);
    await userEvent.click(screen.getByRole('radio', { name: 'Choose sequences…' }));
    // Starts on the primary chromosomes.
    expect(screen.getByRole('checkbox', { name: /^chr1 / })).toBeChecked();
    expect(screen.getByRole('checkbox', { name: /^chrEBV / })).not.toBeChecked();
    expect(screen.getByText(/2 sequences/)).toBeInTheDocument();
    await userEvent.click(screen.getByRole('checkbox', { name: /^chr2 / }));
    expect(screen.getByText(/1 sequence,/)).toBeInTheDocument();
    expect(screen.getByText(/--remote \S+ -n chr1 -o hg38\.fa\.gz$/)).toBeInTheDocument();
    await userEvent.click(screen.getByRole('tab', { name: 'Python' }));
    expect(screen.getByText(/export_fasta\(.*"hg38\.fa\.gz", \["chr1"\]\)/)).toBeInTheDocument();

    await userEvent.click(screen.getByRole('button', { name: /^Download \.fa\.gz/ }));
    await screen.findByText('Saved');
    const exported = vi.mocked(fastaChunks).mock.lastCall?.[0].collection;
    expect(exported?.seqs.map((s) => s.name)).toEqual(['chr1']);
  });

  it('disables sequences it cannot rebuild, and downloads the rest', async () => {
    serveSequences([...PICK_SEQS.slice(0, 2), row('prot1', 'protein')]);
    renderPanel();
    await screen.findByText(/3 sequences/);
    await userEvent.click(screen.getByRole('radio', { name: 'Choose sequences…' }));
    const prot = screen.getByRole('checkbox', { name: /^prot1 / });
    expect(prot).toBeDisabled();
    expect(prot).not.toBeChecked();
    expect(screen.getByText(/protein alphabet, not supported here/)).toBeInTheDocument();
    expect(screen.queryByRole('alert')).toBeNull();
    await userEvent.click(screen.getByRole('button', { name: 'Select all' }));
    expect(screen.getByText(/2 sequences/)).toBeInTheDocument();
    await userEvent.click(screen.getByRole('button', { name: /^Download \.fa\.gz/ }));
    await screen.findByText('Saved');
    const exported = vi.mocked(fastaChunks).mock.lastCall?.[0].collection;
    expect(exported?.seqs.map((s) => s.name)).toEqual(['chr1', 'chr2']);
  });

  it('explains, without saving, when all sequences include one it cannot rebuild', async () => {
    serveSequences([row('prot1', 'protein')]);
    renderPanel();
    await screen.findByText(/1 sequence,/);
    await userEvent.click(screen.getByRole('button', { name: /^Download \.fa\.gz/ }));
    expect(await screen.findByRole('alert')).toHaveTextContent(/prot1.*protein alphabet/);
    expect(saveStream).not.toHaveBeenCalled();
  });

  it('stops a download in flight when the page goes away', async () => {
    let signal: AbortSignal | undefined;
    vi.mocked(saveStream).mockImplementationOnce(async (opts) => {
      signal = opts.signal;
      await new Promise(() => {});
    });
    const view = renderPanel();
    await screen.findByText(/\(87 B uncompressed\)/);
    await userEvent.click(screen.getByRole('button', { name: /^Download \.fa\.gz/ }));
    expect(signal?.aborted).toBe(false);
    view.unmount();
    expect(signal?.aborted).toBe(true);
  });

  it('keeps the choice when the tab is left and opened again', async () => {
    renderPanel();
    await screen.findByText(/\(87 B uncompressed\)/);
    await userEvent.click(screen.getByRole('radio', { name: 'Choose sequences…' }));
    await userEvent.click(screen.getByRole('button', { name: 'Toggle' }));
    expect(screen.queryByRole('radio')).toBeNull();
    await userEvent.click(screen.getByRole('button', { name: 'Toggle' }));
    expect(screen.getByRole('radio', { name: 'Choose sequences…' })).toBeChecked();
  });

  it("offers refgenie with this server's install line", async () => {
    renderPanel();
    await screen.findByText(/\(87 B uncompressed\)/);
    await userEvent.click(screen.getByRole('tab', { name: 'refgenie CLI' }));
    expect(screen.getByText(/refgenie pull hg38\/fasta/)).toBeInTheDocument();
    // The same install line as the landing page, for the server's version.
    expect(screen.getByText('pip install "refgenie>=1.0.0a3"')).toBeInTheDocument();
  });
});
