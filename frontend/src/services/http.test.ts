import { describe, expect, it, vi } from 'vitest';
import { ApiClient, ApiError, buildQuery } from './http';

function jsonResponse(body: unknown, status = 200): Response {
  return new Response(JSON.stringify(body), {
    status,
    headers: { 'content-type': 'application/json' },
  });
}

describe('buildQuery', () => {
  it('drops undefined, null and empty string', () => {
    expect(buildQuery({ a: undefined, b: null, c: '', d: 'keep' })).toBe('?d=keep');
  });

  it('comma-joins arrays (the search_fields contract)', () => {
    expect(buildQuery({ search_fields: ['name', 'digest'] })).toBe(
      '?search_fields=name%2Cdigest',
    );
  });

  it('drops an empty array entirely', () => {
    expect(buildQuery({ search_fields: [] })).toBe('');
  });

  it('preserves booleans and zero', () => {
    expect(buildQuery({ expand: true, offset: 0 })).toBe('?expand=true&offset=0');
  });

  it('returns an empty string for no params', () => {
    expect(buildQuery(undefined)).toBe('');
    expect(buildQuery({})).toBe('');
  });
});

describe('ApiClient.url', () => {
  it('normalizes trailing slashes on the base URL', () => {
    const client = new ApiClient({ baseUrl: 'https://api.refgenie.org/v4///' });
    expect(client.url('/genomes')).toBe('https://api.refgenie.org/v4/genomes');
  });

  it('accepts a relative base and a path without a leading slash', () => {
    const client = new ApiClient({ baseUrl: '/v4' });
    expect(client.url('genomes')).toBe('/v4/genomes');
  });

  it('appends the query string, for download anchors', () => {
    const client = new ApiClient({ baseUrl: '/v4' });
    expect(client.url('/relationships/abc', { expand: true })).toBe(
      '/v4/relationships/abc?expand=true',
    );
  });
});

describe('ApiError extraction', () => {
  it('reads a string detail', async () => {
    const fetchImpl = vi.fn(async () => jsonResponse({ detail: 'Genome x not found' }, 404));
    const client = new ApiClient({ baseUrl: '/v4', fetchImpl: fetchImpl as never });
    const error = (await client.get('/genomes/x').catch((e) => e)) as ApiError;
    expect(error).toBeInstanceOf(ApiError);
    expect(error.status).toBe(404);
    expect(error.isNotFound).toBe(true);
    expect(error.detail).toBe('Genome x not found');
  });

  it('flattens a validation-error array detail', async () => {
    const fetchImpl = vi.fn(async () =>
      jsonResponse(
        { detail: [{ loc: ['query', 'limit'], msg: 'must be <= 1000', type: 'value_error' }] },
        422,
      ),
    );
    const client = new ApiClient({ baseUrl: '/v4', fetchImpl: fetchImpl as never });
    const error = (await client.get('/genomes').catch((e) => e)) as ApiError;
    expect(error.status).toBe(422);
    expect(error.detail).toBe('query.limit: must be <= 1000');
  });

  it('truncates an HTML error body', async () => {
    const body = `<html><body>${'x'.repeat(500)}</body></html>`;
    const fetchImpl = vi.fn(
      async () => new Response(body, { status: 502, statusText: 'Bad Gateway' }),
    );
    const client = new ApiClient({ baseUrl: '/v4', fetchImpl: fetchImpl as never });
    const error = (await client.get('/genomes').catch((e) => e)) as ApiError;
    expect(error.detail.length).toBeLessThanOrEqual(201);
    expect(error.detail.startsWith('<html>')).toBe(true);
  });

  it('falls back to statusText on an empty body', async () => {
    const fetchImpl = vi.fn(
      async () => new Response('', { status: 500, statusText: 'Internal Server Error' }),
    );
    const client = new ApiClient({ baseUrl: '/v4', fetchImpl: fetchImpl as never });
    const error = (await client.get('/genomes').catch((e) => e)) as ApiError;
    expect(error.detail).toBe('Internal Server Error');
  });

  it('reports a transport failure as status 0', async () => {
    const fetchImpl = vi.fn(async () => {
      throw new TypeError('Failed to fetch');
    });
    const client = new ApiClient({ baseUrl: '/v4', fetchImpl: fetchImpl as never });
    const error = (await client.get('/genomes').catch((e) => e)) as ApiError;
    expect(error.status).toBe(0);
    expect(error.isNetworkError).toBe(true);
  });
});

describe('AbortSignal propagation', () => {
  it('aborts the underlying request when the caller signal aborts', async () => {
    const controller = new AbortController();
    const fetchImpl = vi.fn((_url: string, init?: RequestInit) => {
      return new Promise<Response>((_resolve, reject) => {
        init?.signal?.addEventListener('abort', () =>
          reject(new DOMException('Aborted', 'AbortError')),
        );
      });
    });
    const client = new ApiClient({ baseUrl: '/v4', fetchImpl: fetchImpl as never });
    const promise = client.get('/genomes', undefined, { signal: controller.signal });
    controller.abort();
    await expect(promise).rejects.toThrow();
    expect(fetchImpl).toHaveBeenCalledTimes(1);
  });
});

describe('error envelope', () => {
  it('reads code, message and field from the local envelope', async () => {
    const fetchImpl = vi.fn(async () =>
      jsonResponse(
        {
          ok: false,
          error: {
            code: 'missing_build_input',
            message: 'fasta is required',
            detail: 'no such file',
            field: 'files.fasta',
          },
        },
        400,
      ),
    );
    const client = new ApiClient({ baseUrl: '/v1', fetchImpl: fetchImpl as never });
    const error = (await client
      .mutate('/actions/build', { method: 'POST', action: 'build', body: {} })
      .catch((e) => e)) as ApiError;

    expect(error.code).toBe('missing_build_input');
    expect(error.field).toBe('files.fasta');
    expect(error.detail).toBe('fasta is required');
  });

  it('reports an unknown code when the server sent no envelope', async () => {
    const fetchImpl = vi.fn(async () => jsonResponse({ detail: 'nope' }, 404));
    const client = new ApiClient({ baseUrl: '/v1', fetchImpl: fetchImpl as never });
    const error = (await client.get('/jobs/x').catch((e) => e)) as ApiError;
    expect(error.code).toBe('unknown');
    expect(error.field).toBeNull();
  });
});

describe('ApiClient.mutate', () => {
  it('always sets the action header, with and without a body', async () => {
    const fetchImpl = vi.fn(async (_url: string, _init?: RequestInit) =>
      jsonResponse({ ok: true }),
    );
    const client = new ApiClient({ baseUrl: '/v1', fetchImpl: fetchImpl as never });

    await client.mutate('/actions/pull', { method: 'POST', action: 'pull', body: { a: 1 } });
    await client.mutate('/actions/assets/x', { method: 'DELETE', action: 'asset.delete' });

    for (const [, init] of fetchImpl.mock.calls) {
      const headers = init?.headers as Record<string, string>;
      expect(headers['X-Refgenie-Action']).toBe('1');
    }
    // The body carries a JSON content type; the header-only DELETE does not.
    const first = fetchImpl.mock.calls[0][1];
    expect((first?.headers as Record<string, string>)['content-type']).toBe(
      'application/json',
    );
  });
});
