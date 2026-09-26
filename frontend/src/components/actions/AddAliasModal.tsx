/**
 * Point a new name at a genome.
 *
 * Shared by the Manage page's alias panel and the Aliases page, so adding an
 * alias is offered wherever aliases are listed rather than only on Manage.
 * The backend action (`alias.set`) is the same either way.
 */

import { useEffect, useState } from 'react';
import { useLocalApiClient } from '../../hooks/useApiClient';
import { useToast } from '../../hooks/useToast';
import { setAlias } from '../../services/actions';
import { INVALIDATION_FOR_ACTION, invalidate } from '../../services/invalidation';
import { ApiError } from '../../services/http';
import { BaseModal } from '../common/BaseModal';
import { FormField } from '../common/FormField';
import { GenomeSelect } from '../build/GenomeSelect';

export interface AddAliasModalProps {
  isOpen: boolean;
  onClose: () => void;
}

export function AddAliasModal({ isOpen, onClose }: AddAliasModalProps) {
  const client = useLocalApiClient();
  const toast = useToast();
  const [alias, setAliasName] = useState('');
  const [genome, setGenome] = useState('');
  const [digest, setDigest] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);

  useEffect(() => {
    if (!isOpen) return;
    setAliasName('');
    setGenome('');
    setDigest(null);
    setError(null);
    setBusy(false);
  }, [isOpen]);

  const submit = async (event: React.FormEvent) => {
    event.preventDefault();
    if (!alias.trim() || !digest) return;
    setBusy(true);
    setError(null);
    try {
      await setAlias(client, { alias: alias.trim(), genome_digest: digest });
      invalidate(INVALIDATION_FOR_ACTION['alias.set']);
      toast.success(`Alias ${alias.trim()} added.`);
      onClose();
    } catch (err) {
      setError(err instanceof ApiError ? err.detail : String(err));
    } finally {
      setBusy(false);
    }
  };

  return (
    <BaseModal isOpen={isOpen} onClose={onClose} title="Add alias" size="md">
      <form className="flex flex-col gap-4" onSubmit={(event) => void submit(event)}>
        <FormField htmlFor="alias-name" label="Alias" required>
          <input
            id="alias-name"
            className="rg-field__input"
            type="text"
            value={alias}
            onChange={(event) => setAliasName(event.target.value)}
          />
        </FormField>

        <GenomeSelect
          value={genome}
          onChange={(name, resolved) => {
            setGenome(name);
            setDigest(resolved);
          }}
          error={genome && !digest ? 'Pick a genome from the list.' : undefined}
        />

        {error && (
          <p className="rg-field__error" role="alert">
            {error}
          </p>
        )}

        <BaseModal.Footer>
          <button type="button" className="rg-btn" onClick={onClose}>
            Cancel
          </button>
          <button
            type="submit"
            className="rg-btn rg-btn--primary"
            disabled={busy || !alias.trim() || !digest}
          >
            {busy ? 'Adding…' : 'Add alias'}
          </button>
        </BaseModal.Footer>
      </form>
    </BaseModal>
  );
}
