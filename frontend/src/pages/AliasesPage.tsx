import { useMemo, useState } from 'react';
import { useAliases } from '../hooks/queries/useAliases';
import { useGenomeIndex } from '../hooks/queries/useGenomes';
import { useCapability } from '../hooks/useCapability';
import { AddAliasModal } from '../components/actions/AddAliasModal';
import { DigestChip } from '../components/common/DigestChip';
import { ScopedLink } from '../components/common/ScopedLink';
import { SearchListPage } from '../components/common/SearchListPage';
import { ALIAS_SEARCH_FIELDS } from '../services/resources/aliases';
import type { Column } from '../components/common/DataTable';
import type { AliasPublic } from '../types/api';
import type { GenomeIndex } from '../utils/genomes';

function aliasColumns(index: GenomeIndex | undefined): Array<Column<AliasPublic>> {
  const genome = (alias: AliasPublic) => index?.get(alias.genome_digest);
  const na = <span className="rg-muted">NA</span>;
  return [
    {
      key: 'name',
      header: 'Alias',
      render: (alias) => (
        <ScopedLink
          className="rg-link font-medium"
          to={`/genomes/${alias.genome_digest}`}
        >
          {alias.name}
        </ScopedLink>
      ),
    },
    {
      key: 'species',
      header: 'Species',
      render: (alias) => {
        const g = genome(alias);
        return g?.species_name ?? g?.common_name ?? na;
      },
    },
    {
      key: 'assembly',
      header: 'Assembly',
      render: (alias) => {
        const g = genome(alias);
        const parts = [g?.assembly_source, g?.assembly_accession].filter(Boolean);
        return parts.length ? parts.join(' · ') : na;
      },
    },
    {
      key: 'assets',
      header: 'Assets',
      align: 'right',
      render: (alias) => {
        const count = genome(alias)?.asset_count ?? 0;
        return count > 0 ? (
          <span className="font-semibold">{count}</span>
        ) : (
          <span className="rg-muted">0</span>
        );
      },
    },
    {
      key: 'digest',
      header: 'Genome digest',
      render: (alias) => <DigestChip digest={alias.genome_digest} />,
    },
  ];
}

export function AliasesPage() {
  const canWrite = useCapability('aliases_write');
  const [addOpen, setAddOpen] = useState(false);
  const genomes = useGenomeIndex();
  const columns = useMemo(() => aliasColumns(genomes.data), [genomes.data]);

  return (
    <SearchListPage
      title="Aliases"
      lede={
        <>
          An alias is a human-readable name, such as{' '}
          <code className="rg-code rg-code--inline">hg38</code> or{' '}
          <code className="rg-code rg-code--inline">mm10</code>, that resolves to a
          genome digest. You type the alias; refgenie stores and compares everything by
          the digest, and one genome can answer to several aliases.
        </>
      }
      actions={
        canWrite && (
          <button
            type="button"
            className="rg-btn rg-btn--primary"
            onClick={() => setAddOpen(true)}
          >
            Add alias
          </button>
        )
      }
      useList={useAliases}
      searchFields={ALIAS_SEARCH_FIELDS}
      columns={columns}
      compact
      rowKey={(alias) => alias.name}
      emptyMessage="No aliases defined."
    >
      {canWrite && <AddAliasModal isOpen={addOpen} onClose={() => setAddOpen(false)} />}
    </SearchListPage>
  );
}
