import { describe, expect, it } from 'vitest';
import { screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { RemotePage } from './RemotePage';
import { renderWithProviders, serverConfig } from '../test/renderWithProviders';

describe('RemotePage', () => {
  it('renders NotAvailablePage when remote_browse is off', () => {
    renderWithProviders(<RemotePage />, { config: serverConfig });
    expect(screen.getByText('Remote browse')).toBeInTheDocument();
    expect(screen.getByText(/remote_browse/)).toBeInTheDocument();
  });

  it('lists servers with their reachability and error text', async () => {
    renderWithProviders(<RemotePage />);
    expect(await screen.findByText('https://api.refgenie.org')).toBeInTheDocument();
    // Never a silent empty list: the unreachable server shows its error inline.
    expect(screen.getByText('Connection refused')).toBeInTheDocument();
  });

  it('marks a remote genome that is already local', async () => {
    renderWithProviders(<RemotePage />);
    expect(await screen.findByRole('button', { name: 'hg38' })).toBeInTheDocument();
    expect(screen.getAllByText('local').length).toBeGreaterThan(0);
    expect(screen.getAllByText('not local').length).toBeGreaterThan(0);
  });

  it('loads remote assets when a genome is expanded', async () => {
    renderWithProviders(<RemotePage />);
    await userEvent.click(await screen.findByRole('button', { name: 'hg38' }));
    expect(await screen.findByText('fasta:default')).toBeInTheDocument();
  });
});
