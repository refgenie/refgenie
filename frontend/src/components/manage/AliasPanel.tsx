/**
 * Alias curation.
 *
 * The removal path is the interesting one: `alias.remove` has no last-alias
 * guard, so removing a genome's only alias leaves it reachable by digest alone.
 * The confirm still proceeds — the user is simply told the truth first.
 */

import { useMemo, useState } from 'react';
import { Link } from 'react-router-dom';
import { useAliases } from '../../hooks/queries/useAliases';
import { useLocalApiClient } from '../../hooks/useApiClient';
import { useCapability } from '../../hooks/useCapability';
import { useToast } from '../../hooks/useToast';
import { removeAlias, setAlias } from '../../services/actions';
import { INVALIDATION_FOR_ACTION, invalidate } from '../../services/invalidation';
import { ApiError } from '../../services/http';
import { DataTable } from '../common/DataTable';
import { DigestChip } from '../common/DigestChip';
import { Pagination } from '../common/Pagination';
import { BaseModal } from '../common/BaseModal';
import { ConfirmModal } from '../common/ConfirmModal';
import { FormField } from '../common/FormField';
import { GenomeSelect } from '../build/GenomeSelect';
import type { Column } from '../common/DataTable';
import type { AliasPublic } from '../../types/api';

const PAGE_SIZE = 25;
/** Enough to count aliases per genome on a local install without an N+1. */
const COUNT_LIMIT = 1000;

export function AliasPanel() {
  const canWrite = useCapability('aliases_write');
  const client = useLocalApiClient();
  const toast = useToast();
  const [offset, setOffset] = useState(0);
  const [addOpen, setAddOpen] = useState(false);
  const [pendingRemoval, setPendingRemoval] = useState<AliasPublic | null>(null);

  const page = useAliases({ offset, limit: PAGE_SIZE });
  const all = useAliases({ limit: COUNT_LIMIT });

  const countByGenome = useMemo(() => {
    const counts = new Map<string, number>();
    for (const alias of all.data?.items ?? []) {
      counts.set(alias.genome_digest, (counts.get(alias.genome_digest) ?? 0) + 1);
    }
    return counts;
  }, [all.data]);

  const [newAlias, setNewAlias] = useState('');
  const [newGenome, setNewGenome] = useState('');
  const [newDigest, setNewDigest] = useState<string | null>(null);
  const [addError, setAddError] = useState<string | null>(null);
  const [adding, setAdding] = useState(false);

  const resetAddForm = () => {
    setNewAlias('');
    setNewGenome('');
    setNewDigest(null);
    setAddError(null);
  };

  const submitAdd = async (event: React.FormEvent) => {
    event.preventDefault();
    if (!newAlias.trim() || !newDigest) return;
    setAdding(true);
    setAddError(null);
    try {
      await setAlias(client, { alias: newAlias.trim(), genome_digest: newDigest });
      invalidate(INVALIDATION_FOR_ACTION['alias.set']);
      toast.success(`Alias ${newAlias.trim()} added.`);
      setAddOpen(false);
      resetAddForm();
    } catch (error) {
      setAddError(error instanceof ApiError ? error.detail : String(error));
    } finally {
      setAdding(false);
    }
  };

  const columns: Array<Column<AliasPublic>> = [
    {
      key: 'name',
      header: 'Alias',
      render: (alias) => (
        <Link className="rg-link font-medium" to={`/genomes/${alias.genome_digest}`}>
          {alias.name}
        </Link>
      ),
    },
    {
      key: 'digest',
      header: 'Genome',
      render: (alias) => <DigestChip digest={alias.genome_digest} />,
    },
    {
      key: 'siblings',
      header: 'Other aliases',
      align: 'right',
      render: (alias) => Math.max(0, (countByGenome.get(alias.genome_digest) ?? 1) - 1),
    },
    ...(canWrite
      ? [
          {
            key: 'actions',
            header: '',
            render: (alias: AliasPublic) => (
              <button
                type="button"
                className="rg-btn rg-btn--sm rg-btn--danger"
                onClick={() => setPendingRemoval(alias)}
              >
                Remove
              </button>
            ),
          },
        ]
      : []),
  ];

  const removalIsLast =
    pendingRemoval !== null && (countByGenome.get(pendingRemoval.genome_digest) ?? 1) <= 1;

  return (
    <section className="flex flex-col gap-4">
      <header className="flex items-center justify-between gap-4 flex-wrap">
        <h2 className="text-xl font-semibold">Aliases</h2>
        {canWrite && (
          <button
            type="button"
            className="rg-btn rg-btn--primary rg-btn--sm"
            onClick={() => {
              resetAddForm();
              setAddOpen(true);
            }}
          >
            Add alias
          </button>
        )}
      </header>

      <DataTable
        caption="Aliases"
        columns={columns}
        rows={page.data?.items}
        rowKey={(alias) => alias.name}
        loading={page.isPending}
        error={page.error}
        onRetry={() => page.refetch()}
      />

      <Pagination pagination={page.data?.pagination} onOffsetChange={setOffset} />

      <BaseModal
        isOpen={addOpen}
        onClose={() => setAddOpen(false)}
        title="Add alias"
        size="md"
      >
        <form className="flex flex-col gap-4" onSubmit={(event) => void submitAdd(event)}>
          <FormField htmlFor="alias-name" label="Alias" required>
            <input
              id="alias-name"
              className="rg-field__input"
              type="text"
              value={newAlias}
              onChange={(event) => setNewAlias(event.target.value)}
            />
          </FormField>

          <GenomeSelect
            value={newGenome}
            onChange={(alias, digest) => {
              setNewGenome(alias);
              setNewDigest(digest);
            }}
            error={newGenome && !newDigest ? 'Pick a genome from the list.' : undefined}
          />

          {addError && (
            <p className="rg-field__error" role="alert">
              {addError}
            </p>
          )}

          <BaseModal.Footer>
            <button type="button" className="rg-btn" onClick={() => setAddOpen(false)}>
              Cancel
            </button>
            <button
              type="submit"
              className="rg-btn rg-btn--primary"
              disabled={adding || !newAlias.trim() || !newDigest}
            >
              {adding ? 'Adding…' : 'Add alias'}
            </button>
          </BaseModal.Footer>
        </form>
      </BaseModal>

      <ConfirmModal
        isOpen={pendingRemoval !== null}
        onClose={() => setPendingRemoval(null)}
        title="Remove alias"
        destructive
        confirmLabel="Remove alias"
        body={
          <>
            <p className="mb-2">
              Remove <code className="rg-code rg-code--inline">{pendingRemoval?.name}</code>?
            </p>
            {removalIsLast && (
              <p className="rg-banner rg-banner--warning">
                This is the genome&rsquo;s last alias. Removing it leaves the genome reachable
                only by digest.
              </p>
            )}
          </>
        }
        onConfirm={async () => {
          if (!pendingRemoval) return;
          await removeAlias(client, pendingRemoval.name);
          invalidate(INVALIDATION_FOR_ACTION['alias.remove']);
          toast.success(`Alias ${pendingRemoval.name} removed.`);
          setPendingRemoval(null);
        }}
      />
    </section>
  );
}
