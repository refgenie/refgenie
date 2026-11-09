import { describe, expect, it, vi } from 'vitest';
import { act, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { HttpResponse, http } from 'msw';
import { CompareFasta } from './CompareFasta';
import { server } from '../../test/server';
import { API } from '../../test/handlers';
import { renderWithProviders } from '../../test/renderWithProviders';
import type { DigestWorkerReply, WorkerLike } from '../../vendor/fastaDigest/fastaDigestClient.ts';
import type { FastaDigestResult } from '../../vendor/fastaDigest/fastaDigest.ts';
import type { SequenceRow } from './types';

const DIGEST = 'GenomeDigestGenomeDigestGenomeDi';

function seq(name: string, length: number, sha512t24u: string): SequenceRow {
  return {
    name,
    length,
    sha512t24u,
    md5: null,
    alphabet: 'dna2bit',
    description: null,
    aliases: [],
  };
}

const SEQUENCES = [
  seq('chr1', 100, 's1'),
  seq('chr2', 80, 's2'),
  seq('chrM', 16, 'sM'),
];

function result(digest: string, rows: SequenceRow[]): FastaDigestResult {
  return {
    digest,
    n_sequences: rows.length,
    names_digest: 'n',
    sequences_digest: 's',
    lengths_digest: 'l',
    sequences: rows.map((r) => ({
      name: r.name,
      length: r.length,
      sha512t24u: r.sha512t24u,
      md5: '',
      alphabet: 'dna2bit',
      description: null,
    })),
  };
}

const STATS = { chunks: 1, totalBytes: 1000, elapsedMs: 1, avgChunkSize: 1000 };

/** A fake worker the test drives by hand. */
function fakeWorker() {
  const worker: WorkerLike & {
    files: File[];
    send: (reply: DigestWorkerReply) => void;
  } = {
    files: [],
    onmessage: null,
    onerror: null,
    terminate: vi.fn(),
    postMessage(message) {
      worker.files.push(message.file);
    },
    send(reply) {
      act(() => worker.onmessage?.({ data: reply } as MessageEvent<DigestWorkerReply>));
    },
  };
  return worker;
}

function setup() {
  const worker = fakeWorker();
  renderWithProviders(
    <CompareFasta digest={DIGEST} sequences={SEQUENCES} workerFactory={() => worker} />,
  );
  return worker;
}

const fasta = (name = 'mine.fa.gz', bytes = 1000) =>
  new File([new Uint8Array(bytes)], name, { type: 'application/gzip' });

describe('CompareFasta', () => {
  it('says the file never leaves the browser', () => {
    setup();
    expect(
      screen.getByRole('heading', { name: 'Compare your own FASTA' }),
    ).toBeInTheDocument();
    expect(screen.getByText(/never leaves your computer/)).toBeInTheDocument();
  });

  it('digests a picked file with progress, then shows an identical match', async () => {
    const worker = setup();
    await userEvent.upload(screen.getByLabelText('FASTA file to compare'), fasta());
    expect(worker.files[0].name).toBe('mine.fa.gz');
    worker.send({ type: 'progress', bytesRead: 500, totalBytes: 1000 });
    expect(
      screen.getByRole('progressbar', { name: /mine\.fa\.gz: 50%/ }),
    ).toBeInTheDocument();
    worker.send({ type: 'result', result: result(DIGEST, SEQUENCES), stats: STATS });
    expect(
      await screen.findByRole('heading', { name: 'Your file is this genome' }),
    ).toBeInTheDocument();
    expect(screen.getByText('Same coordinates')).toBeInTheDocument();
    expect(screen.getByTitle(DIGEST)).toBeInTheDocument();
    expect(
      screen.getByRole('columnheader', { name: 'In your file' }),
    ).toBeInTheDocument();
    expect(worker.terminate).toHaveBeenCalled();
  });

  it('explains a renamed and a dropped sequence', async () => {
    server.use(
      http.get(`${API}/genomes/:digest`, () =>
        HttpResponse.json({ detail: 'not found' }, { status: 404 }),
      ),
    );
    const worker = setup();
    await userEvent.upload(
      screen.getByLabelText('FASTA file to compare'),
      fasta('mine.fa'),
    );
    worker.send({
      type: 'result',
      result: result('FileDigest', [seq('chr1', 100, 's1'), seq('MT', 16, 'sM')]),
      stats: STATS,
    });
    expect(
      await screen.findByRole('heading', { name: 'Your file is part of this genome' }),
    ).toBeInTheDocument();
    expect(
      screen.getByText('Your file is missing 1 sequence of this genome (chr2).'),
    ).toBeInTheDocument();
    expect(
      screen.getByText(
        '1 shared sequence has a different name in your file (chrM is MT).',
      ),
    ).toBeInTheDocument();
    expect(
      screen.getByText('The coordinate systems partly overlap.'),
    ).toBeInTheDocument();
    expect(screen.queryByText(/matches/)).not.toBeInTheDocument();
  });

  it('links another genome in refgenie with the same digest', async () => {
    server.use(
      http.get(`${API}/genomes/:digest`, ({ params }) =>
        HttpResponse.json({
          digest: params.digest,
          synonyms: null,
          species_name: null,
        }),
      ),
      http.get(`${API}/aliases`, ({ request }) => {
        const genome = new URL(request.url).searchParams.get('genome_digest');
        return HttpResponse.json({
          results: [
            { name: 'hg38_alt_long', genome_digest: genome },
            { name: 'hg38_alt', genome_digest: genome },
          ],
          pagination: { page: 0, page_size: 100, total: 2 },
        });
      }),
    );
    const worker = setup();
    await userEvent.upload(screen.getByLabelText('FASTA file to compare'), fasta());
    worker.send({
      type: 'result',
      result: result('OtherGenome', [seq('chr1', 100, 's1')]),
      stats: STATS,
    });
    expect(await screen.findByRole('link', { name: 'hg38_alt' })).toHaveAttribute(
      'href',
      '/genomes/OtherGenome',
    );
  });

  it('cancels a running digest', async () => {
    const worker = setup();
    await userEvent.upload(screen.getByLabelText('FASTA file to compare'), fasta());
    await userEvent.click(screen.getByRole('button', { name: 'Cancel' }));
    expect(worker.terminate).toHaveBeenCalled();
    expect(screen.getByText('Drop a FASTA file here')).toBeInTheDocument();
  });

  it('stops a running digest when it unmounts', async () => {
    const worker = fakeWorker();
    const view = renderWithProviders(
      <CompareFasta digest={DIGEST} sequences={SEQUENCES} workerFactory={() => worker} />,
    );
    await userEvent.upload(screen.getByLabelText('FASTA file to compare'), fasta());
    view.unmount();
    expect(worker.terminate).toHaveBeenCalled();
    // A late reply from the stopped worker is ignored.
    worker.send({ type: 'progress', bytesRead: 10, totalBytes: 1000 });
  });

  it('shows a digest error', async () => {
    const worker = setup();
    await userEvent.upload(screen.getByLabelText('FASTA file to compare'), fasta());
    worker.send({
      type: 'error',
      message: 'No sequences found. Is this a FASTA file?',
    });
    expect(await screen.findByRole('alert')).toHaveTextContent(
      'mine.fa.gz: No sequences found. Is this a FASTA file?',
    );
  });

  it('refuses a file that is not FASTA without starting a worker', async () => {
    const worker = setup();
    await userEvent.upload(
      screen.getByLabelText('FASTA file to compare'),
      fasta('reads.bam'),
      {
        applyAccept: false,
      },
    );
    expect(screen.getByRole('alert')).toHaveTextContent(/Choose a FASTA file/);
    expect(worker.files).toHaveLength(0);
  });
});
