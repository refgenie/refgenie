/**
 * Ingest a FASTA into the RefgetStore and (optionally) build its fasta asset.
 *
 * Included because without it a fresh local install has nothing to build
 * *against*: every other flow on this page presupposes a genome. It produces a
 * job like any other long operation, so it reuses the whole console for free.
 */

import { useEffect, useState } from 'react';
import { useLocalApiClient } from '../../hooks/useApiClient';
import { useToast } from '../../hooks/useToast';
import { initGenome } from '../../services/actions';
import { ApiError } from '../../services/http';
import { provisionalJob, useJobStore } from '../../stores/jobStore';
import { BaseModal } from '../common/BaseModal';
import { FormField } from '../common/FormField';

export interface InitGenomeModalProps {
  isOpen: boolean;
  onClose: () => void;
}

export function InitGenomeModal({ isOpen, onClose }: InitGenomeModalProps) {
  const client = useLocalApiClient();
  const toast = useToast();
  const [fasta, setFasta] = useState('');
  const [aliases, setAliases] = useState('');
  const [description, setDescription] = useState('');
  const [species, setSpecies] = useState('');
  const [buildFasta, setBuildFasta] = useState(true);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (!isOpen) return;
    setFasta('');
    setAliases('');
    setDescription('');
    setSpecies('');
    setBuildFasta(true);
    setBusy(false);
    setError(null);
  }, [isOpen]);

  const aliasList = aliases
    .split(',')
    .map((alias) => alias.trim())
    .filter(Boolean);
  const valid = fasta.trim().length > 0 && aliasList.length > 0;

  const submit = async (event: React.FormEvent) => {
    event.preventDefault();
    if (!valid) return;
    setBusy(true);
    setError(null);
    try {
      const ref = await initGenome(client, {
        fasta: fasta.trim(),
        aliases: aliasList,
        description: description || null,
        species: species || null,
        build_fasta_asset: buildFasta,
      });
      const store = useJobStore.getState();
      if (ref.duplicate) {
        store.focusJob(ref.job_id);
        toast.info('That genome is already being initialized.');
      } else {
        store.registerQueued(
          provisionalJob({
            id: ref.job_id,
            kind: ref.kind ?? 'genome_init',
            status: ref.status ?? 'queued',
            created_at: ref.created_at,
            label: `initialize ${aliasList[0]}`,
            target: {
              genome_digest: null,
              genome_name: aliasList[0],
              asset_group_name: 'fasta',
              asset_name: null,
            },
          }),
        );
        toast.success('Genome initialization queued. Progress is in the job console.');
      }
      onClose();
    } catch (caught) {
      setError(caught instanceof ApiError ? caught.detail : String(caught));
    } finally {
      setBusy(false);
    }
  };

  return (
    <BaseModal isOpen={isOpen} onClose={onClose} title="Initialize genome" size="md">
      <form className="flex flex-col gap-4" onSubmit={(event) => void submit(event)}>
        <FormField
          htmlFor="init-fasta"
          label="FASTA path"
          required
          hint="A path on the machine running refgenie, or a URL it can reach."
        >
          <input
            id="init-fasta"
            className="rg-field__input rg-field__input--mono"
            type="text"
            spellCheck={false}
            placeholder="/absolute/path/to/genome.fa.gz"
            value={fasta}
            onChange={(event) => setFasta(event.target.value)}
          />
        </FormField>

        <FormField
          htmlFor="init-aliases"
          label="Aliases"
          required
          hint="Comma separated. The first becomes the genome's primary name."
        >
          <input
            id="init-aliases"
            className="rg-field__input"
            type="text"
            placeholder="hg38, GRCh38"
            value={aliases}
            onChange={(event) => setAliases(event.target.value)}
          />
        </FormField>

        <FormField htmlFor="init-species" label="Species">
          <input
            id="init-species"
            className="rg-field__input"
            type="text"
            value={species}
            onChange={(event) => setSpecies(event.target.value)}
          />
        </FormField>

        <FormField htmlFor="init-description" label="Description">
          <textarea
            id="init-description"
            className="rg-field__input"
            rows={2}
            value={description}
            onChange={(event) => setDescription(event.target.value)}
          />
        </FormField>

        <label className="rg-check" htmlFor="init-build-fasta">
          <input
            id="init-build-fasta"
            type="checkbox"
            checked={buildFasta}
            onChange={(event) => setBuildFasta(event.target.checked)}
          />
          <span>Build the fasta asset as well</span>
        </label>

        {error && (
          <p className="rg-field__error" role="alert">
            {error}
          </p>
        )}

        <BaseModal.Footer>
          <button type="button" className="rg-btn" onClick={onClose} disabled={busy}>
            Cancel
          </button>
          <button type="submit" className="rg-btn rg-btn--primary" disabled={!valid || busy}>
            {busy ? 'Submitting…' : 'Initialize'}
          </button>
        </BaseModal.Footer>
      </form>
    </BaseModal>
  );
}
