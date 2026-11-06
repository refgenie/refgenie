/**
 * The HTTP layer. Clients are built from the base URL resolved at startup, so
 * there is no hard-coded API base anywhere in `src/`.
 */

import type { ActionVerb } from './contracts';

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
   * Stable error code from `refgenie/server/errors.py`, or `'unknown'`. Branch on
   * this, not on the message text.
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
 * Turn a non-2xx response into an ApiError. Reads the local error envelope, then
 * FastAPI's `detail` (string or 422 array), else the truncated raw body.
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
      const plain = body as { detail?: unknown };
      if (envelope && typeof envelope === 'object') {
        code = asString(envelope.code);
        field = asString(envelope.field) ?? null;
        detail = asString(envelope.message) ?? asString(envelope.detail) ?? '';
      } else if (typeof plain?.detail === 'string') {
        detail = plain.detail;
      } else if (Array.isArray(plain?.detail)) {
        detail = formatValidationDetail(plain.detail as ValidationErrorItem[]);
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
 * Serialize query parameters: drop `undefined`, `null` and `''`, and join arrays
 * with `,` (as `server/search.py` expects for `search_fields`).
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
 * Header every state-changing request must carry. A cross-origin form post or
 * `<img src>` cannot set it without a CORS preflight, which the server declines.
 */
export const ACTION_HEADER = 'X-Refgenie-Action';

/** The server checks only that the header is present, not its value. */
export const ACTION_HEADER_VALUE = '1';

export type MutateMethod = 'POST' | 'PUT' | 'PATCH' | 'DELETE';

export interface MutateOptions {
  method: MutateMethod;
  action: ActionVerb;
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

  /** Full URL for a path; public so download links can be plain `<a href>`. */
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

  /** Send a state-changing request with the `X-Refgenie-Action` header. */
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
      // Rethrow caller aborts as-is so the cache can tell cancellation from failure.
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
