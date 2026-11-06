/**
 * Fake `EventSourceLike` for tests. Records its URL, keeps listeners per named
 * event, and lets a test call `open()`, `fail()` and `close()`.
 */

import type { EventSourceLike } from '../hooks/useJobEvents';

export const createdSources: FakeEventSource[] = [];

export class FakeEventSource implements EventSourceLike {
  readonly url: string;
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
  }

  open(): void {
    this.onopen?.(new Event('open'));
  }

  fail(): void {
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
