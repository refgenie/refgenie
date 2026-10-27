import { describe, expect, it } from 'vitest';
import { screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { HttpResponse, http } from 'msw';
import { RemotePage } from './RemotePage';
import { renderWithProviders, serverConfig } from '../test/renderWithProviders';
import { server } from '../test/server';
import { LOCAL_API } from '../test/handlers';
import { useInvalidationBridge } from '../hooks/useInvalidationBridge';

/** What `AppLayout` mounts in the real app, so an action refreshes the page. */
function InvalidationBridge() {
  useInvalidationBridge();
  return null;
}

describe('RemotePage', () => {
  it('renders NotAvailablePage when remote_browse is off', () => {
    renderWithProviders(<RemotePage />, { config: serverConfig });
    expect(screen.getByText('Remote browse')).toBeInTheDocument();
    expect(screen.getByText(/remote_browse/)).toBeInTheDocument();
    // The gated branch returns before the page's MiniHero, which is what owns
    // the title everywhere else, so NotAvailablePage has to set it itself.
    expect(document.title).toBe('Remote browse · refgenie server');
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

  it('subscribes from the empty state and browses the new server', async () => {
    const NEW_SERVER = 'http://127.0.0.1:8123';
    let subscriptions: string[] = [];
    let genomesQueriedFor: string | null = null;
    server.use(
      http.get(`${LOCAL_API}/remote/servers`, () =>
        HttpResponse.json({
          servers: subscriptions.map((url) => ({
            url,
            subscribed: true,
            reachable: true,
            error: null,
          })),
        }),
      ),
      http.post(`${LOCAL_API}/actions/subscriptions`, async ({ request }) => {
        const body = (await request.json()) as { server_urls: string[] };
        subscriptions = [...subscriptions, ...body.server_urls];
        return HttpResponse.json({
          ok: true,
          message: 'Subscribed',
          data: { subscriptions },
        });
      }),
      http.get(`${LOCAL_API}/remote/genomes`, ({ request }) => {
        genomesQueriedFor = new URL(request.url).searchParams.get('server_url');
        return HttpResponse.json([]);
      }),
    );

    renderWithProviders(
      <>
        <InvalidationBridge />
        <RemotePage />
      </>,
    );
    expect(await screen.findByText(/No servers subscribed/)).toBeInTheDocument();

    await userEvent.type(screen.getByLabelText('Subscribe to a server'), NEW_SERVER);
    await userEvent.click(screen.getByRole('button', { name: 'Subscribe' }));

    // The list refetches off the invalidation, with no manual reload, and the
    // server just added is the one the catalogue is now querying.
    expect(await screen.findByRole('button', { name: /127\.0\.0\.1:8123/ })).toBeInTheDocument();
    await waitFor(() => expect(genomesQueriedFor).toBe(NEW_SERVER));
  });

  it('loads remote assets when a genome is expanded', async () => {
    renderWithProviders(<RemotePage />);
    await userEvent.click(await screen.findByRole('button', { name: 'hg38' }));
    expect(await screen.findByText('fasta:default')).toBeInTheDocument();
  });
});
