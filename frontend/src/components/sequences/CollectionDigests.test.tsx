import { describe, expect, it } from 'vitest';
import { screen, within } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { HttpResponse, http } from 'msw';
import { CollectionDigests } from './CollectionDigests';
import { server } from '../../test/server';
import { localConfig, renderWithProviders } from '../../test/renderWithProviders';

const SELF = 'self-digest';
const OTHER = 'other-digest';

const config = {
  ...localConfig,
  seqcol_services: [{ store: 'main', url: '/seqcol/main', refget_store_url: null }],
};

describe('CollectionDigests', () => {
  it('opens a modal listing the genomes that share an attribute', async () => {
    server.use(
      http.get('*/seqcol/main/list/collection', () =>
        HttpResponse.json({ results: [SELF, OTHER] }),
      ),
    );
    renderWithProviders(
      <CollectionDigests digest={SELF} level1={{ names: 'names-digest' }} />,
      { config },
    );
    expect(screen.queryByRole('dialog')).toBeNull();

    await userEvent.click(screen.getByRole('button', { name: 'Genomes that share this' }));
    const dialog = screen.getByRole('dialog', { name: 'Genomes with the same names' });
    expect(await within(dialog).findByText(/1 other genome has the same/)).toBeInTheDocument();
    expect(within(dialog).getByRole('link')).toHaveAttribute(
      'href',
      `/genomes/${OTHER}?tab=details`,
    );

    await userEvent.keyboard('{Escape}');
    expect(screen.queryByRole('dialog')).toBeNull();
  });
});
