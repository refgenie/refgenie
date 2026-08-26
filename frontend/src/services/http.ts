/**
 * The one HTTP layer. There is deliberately no module-level API base constant
 * anywhere in `src/`: every client is constructed from the config resolved at
 * startup, so a second backend is `createClient(otherUrl)` and nothing else
 * changes.
 */

export type QueryValue =
  | string
  | number
  | boolean
  | string[]
  | null
  | undefined;

export type QueryParams = Record<string, QueryValue>;

interface ValidationErrorItem {
  loc?: unknown[];
  msg?: string;
  type?: string;
}

/** The local error envelope: `{"ok": false, "error": {...}}`. */
interface ErrorEnvelope {
  ok?: boolean;
  error?: {
    code?: unknown;
    message?: unknown;
    detail?: unknown;
    field?: unknown;
  };
}

/** Stand-in when the server did not classify the failure. */
export const UNKNOWN_ERROR_CODE = 'unknown';

/** Normalized error for every non-2xx response and every transport failure. */
export class ApiError extends Error {
  /** HTTP status, or 0 when the request never got a response. */
  readonly status: number;
  readonly detail: string;
  readonly url: string;
  readonly body?: unknown;
  /**
   * Stable snake_case classification from `refgenie/server/errors.py`, or
   * `'unknown'` when the response carried no envelope. An error boundary or a
   * presentation map branches on this, never on the message text.
   */
  readonly code: string;
  /** Dotted field path, present only on input-validation failures. */
  readonly field: string | null;

  constructor(args: {
    status: number;
    detail: string;
    url: string;
    body?: unknown;
    code?: string;
    field?: string | null;
  }) {
    super(args.detail);
    this.name = 'ApiError';
    this.status = args.status;
    this.detail = args.detail;
    this.url = args.url;
    this.body = args.body;
    this.code = args.code ?? UNKNOWN_ERROR_CODE;
    this.field = args.field ?? null;
  }

  get isNotFound(): boolean {
    return this.status === 404;
  }

  /** True when the request never reached the server (offline, CORS, DNS). */
  get isNetworkError(): boolean {
    return this.status === 0;
  }
}

function formatValidationDetail(items: ValidationErrorItem[]): string {
  return items
    .map((item) => {
      const loc = Array.isArray(item.loc) ? item.loc.join('.') : '';
      const msg = item.msg ?? 'invalid value';
      return loc ? `${loc}: ${msg}` : msg;
    })
    .join('; ');
}

function truncate(text: string, max = 200): string {
  const trimmed = text.trim();
  return trimmed.length > max ? `${trimmed.slice(0, max)}…` : trimmed;
}

function asString(value: unknown): string | undefined {
  return typeof value === 'string' && value ? value : undefined;
}

/**
 * Turn a non-2xx response into an ApiError.
 *
 * Three body shapes reach this, in precedence order: the local error envelope
 * `{"ok": false, "error": {code, message, detail, field}}`; FastAPI's `detail`
 * (a string, or the 422 validation array); and anything else — a proxy or a
 * static host may return HTML — which is truncated verbatim.
 */
export async function parseError(response: Response, url?: string): Promise<ApiError> {
  const target = url ?? response.url ?? '';
  let raw = '';
  try {
    raw = await response.text();
  } catch {
    raw = '';
  }

  let body: unknown;
  let detail = '';
  let code: string | undefined;
  let field: string | null = null;

  if (raw) {
    try {
      body = JSON.parse(raw);
      const envelope = (body as ErrorEnvelope)?.error;
      const legacy = body as { detail?: unknown };
      if (envelope && typeof envelope === 'object') {
        code = asString(envelope.code);
        field = asString(envelope.field) ?? null;
        detail = asString(envelope.message) ?? asString(envelope.detail) ?? '';
      } else if (typeof legacy?.detail === 'string') {
        detail = legacy.detail;
      } else if (Array.isArray(legacy?.detail)) {
        detail = formatValidationDetail(legacy.detail as ValidationErrorItem[]);
      } else {
        detail = truncate(raw);
      }
    } catch {
      detail = truncate(raw);
    }
  }
  if (!detail) detail = response.statusText || `HTTP ${response.status}`;

  return new ApiError({ status: response.status, detail, url: target, body, code, field });
}

/**
 * Serialize query parameters the way the backend expects.
 *
 * Drops `undefined`, `null` and `''`; joins arrays with `,` (the
 * `search_fields` contract in `utils/search.py::get_search_params`); leaves
 * booleans as `true` / `false`. Callers pass wire names; nothing is renamed.
 */
export function buildQuery(params: QueryParams | undefined): string {
  if (!params) return '';
  const search = new URLSearchParams();
  for (const [key, value] of Object.entries(params)) {
    if (value === undefined || value === null || value === '') continue;
    if (Array.isArray(value)) {
      const joined = value.filter((v) => v !== '').join(',');
      if (joined) search.append(key, joined);
      continue;
    }
    search.append(key, String(value));
  }
  const qs = search.toString();
  return qs ? `?${qs}` : '';
}

export interface ApiClientOptions {
  baseUrl: string;
  fetchImpl?: typeof fetch;
  timeoutMs?: number;
}

export interface RequestInitLite {
  signal?: AbortSignal;
}

/**
 * The custom header every state-changing request must carry. Presence is the
 * whole signal (the server ignores the value); it exists so a cross-origin
 * form post or an `<img src>` cannot trigger a local mutation, because neither
 * can set a custom header without a preflight the server declines.
 */
export const ACTION_HEADER = 'X-Refgenie-Action';

/** Sent as the header value. The server ignores it; only presence matters. */
export const ACTION_HEADER_VALUE = '1';

export type MutateMethod = 'POST' | 'PUT' | 'PATCH' | 'DELETE';

export interface MutateOptions {
  method: MutateMethod;
  /**
   * The verb this call performs (`'pull'`, `'alias.remove'`, …). Required with
   * no default so a caller physically cannot forget that this is a mutation;
   * it is what puts `X-Refgenie-Action` on the wire.
   */
  action: string;
  body?: unknown;
  params?: QueryParams;
  signal?: AbortSignal;
}

const DEFAULT_TIMEOUT_MS = 30_000;

export class ApiClient {
  readonly baseUrl: string;
  private readonly fetchImpl: typeof fetch;
  private readonly timeoutMs: number;

  constructor(opts: ApiClientOptions) {
    this.baseUrl = opts.baseUrl.replace(/\/+$/, '');
    this.fetchImpl = opts.fetchImpl ?? globalThis.fetch.bind(globalThis);
    this.timeoutMs = opts.timeoutMs ?? DEFAULT_TIMEOUT_MS;
  }

  /**
   * Absolute (or root-relative) URL for a path. Download links must be real
   * `<a href>` anchors rather than fetches, so this is public API.
   */
  url(path: string, params?: QueryParams): string {
    const suffix = path.startsWith('/') ? path : `/${path}`;
    return `${this.baseUrl}${suffix}${buildQuery(params)}`;
  }

  get<T>(path: string, params?: QueryParams, init?: RequestInitLite): Promise<T> {
    return this.request<T>('GET', this.url(path, params), undefined, init);
  }

  post<T>(path: string, body?: unknown, init?: RequestInitLite): Promise<T> {
    return this.request<T>('POST', this.url(path), body, init);
  }

  delete<T>(path: string, params?: QueryParams, init?: RequestInitLite): Promise<T> {
    return this.request<T>('DELETE', this.url(path, params), undefined, init);
  }

  /**
   * The mutation half. Every state-changing request goes through here, which
   * is the only place `X-Refgenie-Action` is attached; no component calls
   * `fetch` directly (an ESLint rule enforces that).
   */
  mutate<T>(path: string, opts: MutateOptions): Promise<T> {
    return this.request<T>(
      opts.method,
      this.url(path, opts.params),
      opts.body,
      { signal: opts.signal },
      { [ACTION_HEADER]: ACTION_HEADER_VALUE },
    );
  }

  private async request<T>(
    method: string,
    url: string,
    body?: unknown,
    init?: RequestInitLite,
    extraHeaders?: Record<string, string>,
  ): Promise<T> {
    const controller = new AbortController();
    const timer = setTimeout(() => controller.abort(), this.timeoutMs);
    const onAbort = () => controller.abort();
    init?.signal?.addEventListener('abort', onAbort);

    try {
      const headers: Record<string, string> = { ...extraHeaders };
      if (body !== undefined) headers['content-type'] = 'application/json';
      const response = await this.fetchImpl(url, {
        method,
        signal: controller.signal,
        headers: Object.keys(headers).length ? headers : undefined,
        body: body === undefined ? undefined : JSON.stringify(body),
      });
      if (!response.ok) throw await parseError(response, url);
      if (response.status === 204) return undefined as T;
      const text = await response.text();
      return (text ? JSON.parse(text) : undefined) as T;
    } catch (error) {
      if (error instanceof ApiError) throw error;
      // A caller-initiated abort must propagate as-is so TanStack Query can
      // distinguish cancellation from failure.
      if (init?.signal?.aborted) throw error;
      throw new ApiError({
        status: 0,
        detail: error instanceof Error ? error.message : 'Network request failed',
        url,
      });
    } finally {
      clearTimeout(timer);
      init?.signal?.removeEventListener('abort', onAbort);
    }
  }
}
