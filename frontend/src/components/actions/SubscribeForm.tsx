/**
 * Subscribe to a refgenie server.
 *
 * It lives in `actions/` because two surfaces need it: the Manage page's
 * subscription panel, and the Remote page's server picker — where a user with
 * no subscriptions actually lands, and where a form beats a sentence naming a
 * CLI command.
 *
 * `subscribe` writes the config, it does not probe: an unreachable URL still
 * subscribes, and the reachability badge on the list is what tells the truth
 * afterwards. A duplicate is a silent set-union server-side, so it is caught
 * here and reported on the field rather than looking like it did something.
 *
 * Self-gating on `subscriptions`, like every other action affordance. It
 * invalidates `servers` AND `remote`, so the catalogue picks the new server up
 * without a reload.
 */

import { useId, useState } from 'react';
import { useLocalApiClient } from '../../hooks/useApiClient';
import { useCapability } from '../../hooks/useCapability';
import { useToast } from '../../hooks/useToast';
import { subscribe } from '../../services/actions';
import { INVALIDATION_FOR_ACTION, invalidate } from '../../services/invalidation';
import { ApiError } from '../../services/http';
import { FormField } from '../common/FormField';

export interface SubscribeFormProps {
  /** The current subscriptions, so a duplicate costs no round trip. */
  knownUrls?: readonly string[];
  /** Manage offers "replace existing subscriptions"; the picker does not. */
  allowReset?: boolean;
  label?: string;
  /** The URL that was just subscribed to, for a caller that selects it. */
  onSubscribed?: (url: string) => void;
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
}: SubscribeFormProps) {
  const canSubscribe = useCapability('subscriptions');
  const client = useLocalApiClient();
  const toast = useToast();
  const inputId = useId();
  const resetId = useId();

  const [url, setUrl] = useState('');
  const [reset, setReset] = useState(false);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);

  if (!canSubscribe) return null;

  const submit = async (event: React.FormEvent) => {
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
    setBusy(true);
    setError(null);
    try {
      // Always the list form, and never empty: the model sets min_length=1.
      await subscribe(client, { server_urls: [target], reset });
      invalidate(INVALIDATION_FOR_ACTION.subscribe);
      toast.success(`Subscribed to ${target}.`);
      setUrl('');
      setReset(false);
      onSubscribed?.(target);
    } catch (caught) {
      setError(caught instanceof ApiError ? caught.detail : String(caught));
    } finally {
      setBusy(false);
    }
  };

  return (
    // `noValidate`: the check in `submit` is the one that speaks, and it
    // speaks on the field. A native validation bubble would say something
    // vaguer, in a tooltip that disappears.
    <form className="flex flex-col gap-3" noValidate onSubmit={(event) => void submit(event)}>
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
