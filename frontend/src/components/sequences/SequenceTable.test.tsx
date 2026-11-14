import { describe, expect, it } from 'vitest';
import { render, screen, within } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { SequenceTable } from './SequenceTable';
import { genomeSequencesFixture } from '../../test/fixtures';
import type { SequenceRow } from './types';

const rows = genomeSequencesFixture.sequences;

function bodyRows() {
  const [, body] = screen.getAllByRole('rowgroup');
  return within(body).getAllByRole('row');
}

function makeRows(n: number): SequenceRow[] {
  return Array.from({ length: n }, (_, i) => ({
    name: `chr${i + 1}`,
    length: 1000 + i,
    sha512t24u: `digest${String(i).padStart(26, '0')}`,
    md5: null,
    alphabet: 'dna2bit',
    description: null,
    aliases: [],
  }));
}

describe('SequenceTable', () => {
  it('renders names, formatted lengths and the preferred alias', () => {
    render(<SequenceTable sequences={rows} />);
    expect(screen.getByText('rCRSd')).toBeInTheDocument();
    expect(screen.getByText('extra')).toBeInTheDocument();
    expect(screen.getByText((33138).toLocaleString())).toBeInTheDocument();
    expect(screen.getByText('NC_012920.1')).toBeInTheDocument();
  });

  it('filters rows by alias', async () => {
    render(<SequenceTable sequences={rows} />);
    await userEvent.type(screen.getByLabelText('Filter sequences'), 'nc_0129');
    await screen.findByText('NC_012920.1');
    await expect.poll(() => bodyRows().length).toBe(1);
    expect(screen.queryByText('extra')).not.toBeInTheDocument();
  });

  it('sorts by length, descending on the second click', async () => {
    render(<SequenceTable sequences={rows} />);
    const button = screen.getByRole('button', { name: /Length/ });
    await userEvent.click(button);
    expect(button.closest('th')).toHaveAttribute('aria-sort', 'ascending');
    expect(within(bodyRows()[0]).getByText('extra')).toBeInTheDocument();
    await userEvent.click(button);
    expect(button.closest('th')).toHaveAttribute('aria-sort', 'descending');
    expect(within(bodyRows()[0]).getByText('rCRSd')).toBeInTheDocument();
  });

  it('pages at 50 rows', async () => {
    render(<SequenceTable sequences={makeRows(60)} />);
    expect(bodyRows()).toHaveLength(50);
    expect(screen.queryByText('chr51')).not.toBeInTheDocument();
    await userEvent.click(screen.getByRole('button', { name: 'Next page' }));
    expect(screen.getByText('chr51')).toBeInTheDocument();
    expect(bodyRows()).toHaveLength(10);
  });

  it('shows a titled dash, not NA, for a sequence with no alias', () => {
    render(<SequenceTable sequences={rows} />);
    expect(screen.queryByText('NA')).not.toBeInTheDocument();
    const extra = bodyRows().find((row) => within(row).queryByText('extra'))!;
    const dash = within(extra).getByTitle(/^No alias/);
    expect(dash).toHaveTextContent('–');
  });

  it('explains the digest and alphabet columns', () => {
    render(<SequenceTable sequences={rows} />);
    expect(screen.getByRole('columnheader', { name: 'SHA-512/24u' })).toHaveAttribute(
      'title',
      expect.stringMatching(/same digest means the same sequence/),
    );
    expect(screen.getAllByTitle(/^dna3bit: A, C, G, T and N/).length).toBeGreaterThan(0);
  });
});
