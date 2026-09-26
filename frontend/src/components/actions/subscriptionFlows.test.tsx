/**
 * Subscribing and unsubscribing, the two halves of the same operation.
 *
 * The interesting assertions are the ones about NOT sending a request: a
 * malformed URL and a duplicate are both answered on the field, because the
 * config write would either 422 or silently union and look like a success.
 */

import { screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { HttpResponse, http } from 'msw';
import { describe, expect, it } from 'vitest';
import { SubscribeForm } from './SubscribeForm';
import { UnsubscribeButton } from './UnsubscribeButton';
import { renderWithProviders, serverConfig } from '../../test/renderWithProviders';
import { server } from '../../test/server';
import { LOCAL_API } from '../../test/handlers';

const SERVER_URL = 'http://127.0.0.1:8123';

describe('SubscribeForm', () => {
  it('posts the URL with the action header and clears the field', async () => {
    let headerSeen: string | null = null;
    let body: unknown;
    const subscribed: string[] = [];
    server.use(
      http.post(`${LOCAL_API}/actions/subscriptions`, async ({ request }) => {
        headerSeen = request.headers.get('x-refgenie-action');
        body = await request.json();
        return HttpResponse.json({
          ok: true,
          message: `Subscribed to: ${SERVER_URL}`,
          data: { subscriptions: [SERVER_URL] },
        });
      }),
    );

    renderWithProviders(<SubscribeForm onSubscribed={(url) => subscribed.push(url)} />);

    const field = screen.getByLabelText('Subscribe to a server');
    await userEvent.type(field, SERVER_URL);
    await userEvent.click(screen.getByRole('button', { name: 'Subscribe' }));

    await waitFor(() => expect(subscribed).toEqual([SERVER_URL]));
    expect(headerSeen).not.toBeNull();
    // Always the list form: the request model sets min_length=1.
    expect(body).toEqual({ server_urls: [SERVER_URL], reset: false });
    expect(field).toHaveValue('');
  });

  it('rejects a URL that is not http(s) on the field, without a request', async () => {
    let called = false;
    server.use(
      http.post(`${LOCAL_API}/actions/subscriptions`, () => {
        called = true;
        return HttpResponse.json({ ok: true, message: '', data: null });
      }),
    );

    renderWithProviders(<SubscribeForm />);
    await userEvent.type(screen.getByLabelText('Subscribe to a server'), 'refgenomes.databio.org');
    await userEvent.click(screen.getByRole('button', { name: 'Subscribe' }));

    expect(await screen.findByRole('alert')).toHaveTextContent('Enter a full http(s) URL');
    expect(called).toBe(false);
  });

  it('rejects a duplicate on the field, without a request', async () => {
    let called = false;
    server.use(
      http.post(`${LOCAL_API}/actions/subscriptions`, () => {
        called = true;
        return HttpResponse.json({ ok: true, message: '', data: null });
      }),
    );

    renderWithProviders(<SubscribeForm knownUrls={[SERVER_URL]} />);
    await userEvent.type(screen.getByLabelText('Subscribe to a server'), SERVER_URL);
    await userEvent.click(screen.getByRole('button', { name: 'Subscribe' }));

    expect(await screen.findByRole('alert')).toHaveTextContent('Already subscribed');
    expect(called).toBe(false);
  });

  it("shows the server's own message when the action fails", async () => {
    server.use(
      http.post(`${LOCAL_API}/actions/subscriptions`, () =>
        HttpResponse.json(
          { ok: false, error: { code: 'conflict', message: 'Configuration is read-only.' } },
          { status: 409 },
        ),
      ),
    );

    renderWithProviders(<SubscribeForm />);
    await userEvent.type(screen.getByLabelText('Subscribe to a server'), SERVER_URL);
    await userEvent.click(screen.getByRole('button', { name: 'Subscribe' }));

    expect(await screen.findByRole('alert')).toHaveTextContent('Configuration is read-only.');
  });

  it('renders nothing when the instance cannot subscribe', () => {
    renderWithProviders(<SubscribeForm />, { config: serverConfig });
    expect(screen.queryByRole('button', { name: 'Subscribe' })).not.toBeInTheDocument();
  });
});

describe('UnsubscribeButton', () => {
  it('confirms, then sends a DELETE carrying the URL', async () => {
    let method: string | null = null;
    let headerSeen: string | null = null;
    let body: unknown;
    const dropped: string[] = [];
    server.use(
      http.delete(`${LOCAL_API}/actions/subscriptions`, async ({ request }) => {
        method = request.method;
        headerSeen = request.headers.get('x-refgenie-action');
        body = await request.json();
        return HttpResponse.json({
          ok: true,
          message: `Unsubscribed from: ${SERVER_URL}`,
          data: { subscriptions: [] },
        });
      }),
    );

    renderWithProviders(
      <UnsubscribeButton url={SERVER_URL} onUnsubscribed={(url) => dropped.push(url)} />,
    );

    await userEvent.click(screen.getByRole('button', { name: 'Unsubscribe' }));
    expect(screen.getByText(/Assets already pulled from it stay/)).toBeInTheDocument();

    // Two now carry the name: the trigger, and the modal's confirm button.
    const buttons = screen.getAllByRole('button', { name: 'Unsubscribe' });
    await userEvent.click(buttons[buttons.length - 1]);

    await waitFor(() => expect(dropped).toEqual([SERVER_URL]));
    expect(method).toBe('DELETE');
    expect(headerSeen).not.toBeNull();
    expect(body).toEqual({ server_urls: [SERVER_URL] });
  });

  it('renders nothing when the instance cannot subscribe', () => {
    renderWithProviders(<UnsubscribeButton url={SERVER_URL} />, { config: serverConfig });
    expect(screen.queryByRole('button', { name: 'Unsubscribe' })).not.toBeInTheDocument();
  });
});
