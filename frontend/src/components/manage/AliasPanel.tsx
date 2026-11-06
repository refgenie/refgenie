/**
 * Alias curation. Removing a genome's last alias is allowed, but the confirm
 * warns that it will then be reachable by digest only.
 */

import { useMemo, useState } from 'react';
import { Link } from 'react-router-dom';
import { useAliases } from '../../hooks/queries/useAliases';
import { useLocalApiClient } from '../../hooks/useApiClient';
import { useCapability } from '../../hooks/useCapability';
import { useToast } from '../../hooks/useToast';
import { useResetPastEndPage } from '../../hooks/useResetPastEndPage';
import { removeAlias } from '../../services/actions';
import { useInvalidate } from '../../hooks/useResource';
import { staleAfter } from '../../services/queryKeys';
import { DataTable } from '../common/DataTable';
import { DigestChip } from '../common/DigestChip';
import { Pagination } from '../common/Pagination';
import { ConfirmModal } from '../common/ConfirmModal';
import { AddAliasModal } from '../actions/AddAliasModal';
import type { Column } from '../common/DataTable';
import type { AliasPublic } from '../../types/api';

const PAGE_SIZE = 25;
/** Enough to count aliases per genome on a local install without an N+1. */
const COUNT_LIMIT = 1000;

export function AliasPanel() {
  const canWrite = useCapability('aliases_write');
  const client = useLocalApiClient();
  const toast = useToast();
  const invalidate = useInvalidate();
  const [pageIndex, setPageIndex] = useState(0);
  const [addOpen, setAddOpen] = useState(false);
  const [pendingRemoval, setPendingRemoval] = useState<AliasPublic | null>(null);

  const page = useAliases({ page: pageIndex, pageSize: PAGE_SIZE });
  // Removing the only alias on the last page leaves that page past the end.
  useResetPastEndPage(page.error, pageIndex, setPageIndex);
  const all = useAliases({ pageSize: COUNT_LIMIT });

  const countByGenome = useMemo(() => {
    const counts = new Map<string, number>();
    for (const alias of all.data?.results ?? []) {
      counts.set(alias.genome_digest, (counts.get(alias.genome_digest) ?? 0) + 1);
    }
    return counts;
  }, [all.data]);

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
            onClick={() => setAddOpen(true)}
          >
            Add alias
          </button>
        )}
      </header>

      <DataTable
        caption="Aliases"
        columns={columns}
        rows={page.data?.results}
        rowKey={(alias) => alias.name}
        loading={page.isPending}
        error={page.error}
        onRetry={() => page.refetch()}
      />

      <Pagination pagination={page.data?.pagination} onPageChange={setPageIndex} />

      <AddAliasModal isOpen={addOpen} onClose={() => setAddOpen(false)} />

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
          invalidate(staleAfter['alias.remove']);
          toast.success(`Alias ${pendingRemoval.name} removed.`);
          setPendingRemoval(null);
        }}
      />
    </section>
  );
}
