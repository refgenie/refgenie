/**
 * Subscribe to a refgenie server. Used by the Servers page and the Remote
 * page's server picker.
 *
 * The server does not probe the URL, and silently ignores a duplicate, so
 * duplicates are caught here and reported on the field.
 */

import { useId, useState } from 'react';
import { useLocalApiClient } from '../../hooks/useApiClient';
import { useCapability } from '../../hooks/useCapability';
import { useToast } from '../../hooks/useToast';
import { subscribe } from '../../services/actions';
import { useInvalidate } from '../../hooks/useResource';
import { staleAfter } from '../../services/queryKeys';
import { useAction } from '../../hooks/useAction';
import { FormField } from '../common/FormField';

export interface SubscribeFormProps {
  /** The current subscriptions, so a duplicate costs no round trip. */
  knownUrls?: readonly string[];
  /** Manage offers "replace existing subscriptions"; the picker does not. */
  allowReset?: boolean;
  label?: string;
  /** The URL that was just subscribed to, for a caller that selects it. */
  onSubscribed?: (url: string) => void;
  /** Prefill the field (never submitted on its own). */
  initialUrl?: string;
  /** Focus the field on mount: set when the user arrived to subscribe. */
  autoFocus?: boolean;
}

function isValidUrl(value: string): boolean {
  try {
    const parsed = new URL(value);
    return parsed.protocol === 'http:' || parsed.protocol === 'https:';
  } catch {
    return false;
  }
}

export function SubscribeForm({
  knownUrls = [],
  allowReset = false,
  label = 'Subscribe to a server',
  onSubscribed,
  initialUrl = '',
  autoFocus = false,
}: SubscribeFormProps) {
  const canSubscribe = useCapability('subscriptions');
  const client = useLocalApiClient();
  const toast = useToast();
  const invalidate = useInvalidate();
  const inputId = useId();
  const resetId = useId();

  const [url, setUrl] = useState(initialUrl);
  const [reset, setReset] = useState(false);
  const { busy, error, setError, run } = useAction();

  if (!canSubscribe) return null;

  const submit = (event: React.FormEvent) => {
    event.preventDefault();
    const target = url.trim();
    if (!isValidUrl(target)) {
      setError('Enter a full http(s) URL, for example https://refgenomes.databio.org.');
      return;
    }
    // With `reset` the duplicate is the point: it becomes the only subscription.
    if (!reset && knownUrls.includes(target)) {
      setError('Already subscribed to this server.');
      return;
    }
    void run(async () => {
      // Always the list form, and never empty: the model sets min_length=1.
      await subscribe(client, { server_urls: [target], reset });
      invalidate(staleAfter.subscribe);
      toast.success(`Subscribed to ${target}.`);
      setUrl('');
      setReset(false);
      onSubscribed?.(target);
    });
  };

  return (
    // `noValidate`: `submit` validates and shows the error on the field.
    <form className="flex flex-col gap-3" noValidate onSubmit={submit}>
      <FormField
        htmlFor={inputId}
        label={label}
        error={error ?? undefined}
        hint="For example https://refgenomes.databio.org"
      >
        <input
          id={inputId}
          className="rg-field__input rg-field__input--mono"
          type="url"
          value={url}
          placeholder="https://refgenomes.databio.org"
          autoFocus={autoFocus}
          onChange={(event) => {
            setUrl(event.target.value);
            setError(null);
          }}
        />
      </FormField>

      {allowReset && (
        <label className="rg-check" htmlFor={resetId}>
          <input
            id={resetId}
            type="checkbox"
            checked={reset}
            onChange={(event) => setReset(event.target.checked)}
          />
          <span>and replace existing subscriptions</span>
        </label>
      )}

      <div className="rg-form-actions">
        <button type="submit" className="rg-btn rg-btn--primary" disabled={busy || !url.trim()}>
          {busy ? 'Subscribing…' : 'Subscribe'}
        </button>
      </div>
    </form>
  );
}
