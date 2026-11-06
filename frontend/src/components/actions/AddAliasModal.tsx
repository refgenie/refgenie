/** Add an alias to a genome. Used on the Manage and Aliases pages. */

import { useEffect, useState } from 'react';
import { useLocalApiClient } from '../../hooks/useApiClient';
import { useToast } from '../../hooks/useToast';
import { setAlias } from '../../services/actions';
import { useInvalidate } from '../../hooks/useResource';
import { staleAfter } from '../../services/queryKeys';
import { useAction } from '../../hooks/useAction';
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
  const invalidate = useInvalidate();
  const [alias, setAliasName] = useState('');
  const [genome, setGenome] = useState('');
  const [digest, setDigest] = useState<string | null>(null);
  const { busy, error, setError, run } = useAction();

  useEffect(() => {
    if (!isOpen) return;
    setAliasName('');
    setGenome('');
    setDigest(null);
    setError(null);
  }, [isOpen, setError]);

  const submit = (event: React.FormEvent) => {
    event.preventDefault();
    if (!alias.trim() || !digest) return;
    void run(async () => {
      await setAlias(client, { alias: alias.trim(), genome_digest: digest });
      invalidate(staleAfter['alias.set']);
      toast.success(`Alias ${alias.trim()} added.`);
      onClose();
    });
  };

  return (
    <BaseModal isOpen={isOpen} onClose={onClose} title="Add alias" size="md">
      <form className="flex flex-col gap-4" onSubmit={submit}>
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
