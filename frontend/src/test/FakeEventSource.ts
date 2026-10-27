/**
 * A hand-written `EventSourceLike`, so the component tests never need a
 * backend and never need a real network connection.
 *
 * It records the URL it was constructed with (the reconnect assertions read
 * `?since=`), keeps listeners per event NAME (the stream uses named events, not
 * default `message` frames), and exposes `open()`, `fail()` and `close()` so a
 * test drives the connection lifecycle explicitly.
 */

import type { EventSourceLike } from '../hooks/useJobEvents';

export const createdSources: FakeEventSource[] = [];

export class FakeEventSource implements EventSourceLike {
  readonly url: string;
  readyState = 0;
  onopen: ((event: Event) => void) | null = null;
  onerror: ((event: Event) => void) | null = null;
  closed = false;

  private readonly listeners = new Map<string, Array<(event: MessageEvent) => void>>();

  constructor(url: string) {
    this.url = url;
    createdSources.push(this);
  }

  addEventListener(type: string, listener: (event: MessageEvent) => void): void {
    const list = this.listeners.get(type) ?? [];
    list.push(listener);
    this.listeners.set(type, list);
  }

  close(): void {
    this.closed = true;
    this.readyState = 2;
  }

  open(): void {
    this.readyState = 1;
    this.onopen?.(new Event('open'));
  }

  fail(): void {
    this.readyState = 2;
    this.onerror?.(new Event('error'));
  }

  /** Wraps the payload the way the real thing does: a JSON string on `data`. */
  emit(type: string, payload: unknown): void {
    this.emitRaw(type, JSON.stringify(payload));
  }

  /** For the malformed-frame case. */
  emitRaw(type: string, data: string): void {
    for (const listener of this.listeners.get(type) ?? []) {
      listener({ data } as MessageEvent);
    }
  }
}

export function resetFakeEventSources(): void {
  createdSources.length = 0;
}

export function fakeEventSourceFactory(url: string): EventSourceLike {
  return new FakeEventSource(url);
}
