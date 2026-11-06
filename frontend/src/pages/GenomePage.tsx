import { Link, useParams } from 'react-router-dom';
import { useGenome } from '../hooks/queries/useGenomes';
import { useAssets } from '../hooks/queries/useAssets';
import { useAliases } from '../hooks/queries/useAliases';
import { useArchives } from '../hooks/queries/useServerInfo';
import { useCapability } from '../hooks/useCapability';
import { useUiConfig } from '../hooks/useUiConfig';
import { preferredAlias } from '../utils/genomes';
import { formatDigest } from '../utils/format';
import { Breadcrumbs } from '../components/common/Breadcrumbs';
import { MiniHero } from '../components/layout/MiniHero';
import { CapabilityGate } from '../components/common/CapabilityGate';
import { DeleteGenomeButton } from '../components/actions/DeleteGenomeButton';
import { ErrorState, LoadingState } from '../components/common/states';
import { GenomeSummaryCard } from '../components/genomes/GenomeSummaryCard';
import { FhrPanel } from '../components/genomes/FhrPanel';
import { GenomeAssetTable } from '../components/genomes/GenomeAssetTable';

export function GenomePage() {
  const { digest } = useParams<{ digest: string }>();
  const config = useUiConfig();
  const canReadArchives = useCapability('archives');

  const genome = useGenome(digest);
  const aliases = useAliases({ genome_digest: digest, pageSize: 100 });
  const assets = useAssets({ genome_digest: digest, pageSize: 200 }, { enabled: !!digest });
  const archives = useArchives({ genome_digest: digest }, { enabled: canReadArchives && !!digest });

  // No early return: the head renders in every state.
  const record = genome.data;
  const aliasNames = aliases.data?.results.map((alias) => alias.name) ?? [];
  // A new genome may have no alias; fall back to a shortened digest.
  const primaryAlias = preferredAlias(aliasNames);
  const title = record ? (primaryAlias ?? formatDigest(record.digest)) : 'Genome';
  // The build form needs an alias; the digest is only for the job target.
  const buildHref =
    `/build?genome=${encodeURIComponent(aliasNames[0] ?? '')}` +
    `&genome_digest=${encodeURIComponent(record?.digest ?? '')}`;

  return (
    <div className="flex flex-col gap-8">
      <MiniHero
        title={title}
        // `false` while loading keeps the placeholder 'Genome' out of the tab title.
        documentTitle={record ? undefined : false}
        breadcrumbs={
          <Breadcrumbs items={[{ label: 'Genomes', to: '/genomes' }, { label: title }]} />
        }
        actions={
          record && (
            <>
              <CapabilityGate cap="build">
                <Link className="rg-btn rg-btn--primary" to={buildHref}>
                  Build asset
                </Link>
              </CapabilityGate>
              <DeleteGenomeButton
                digest={record.digest}
                primaryAlias={primaryAlias ?? ''}
                aliasCount={aliasNames.length}
                assetCount={assets.data?.results.length ?? 0}
              />
            </>
          )
        }
      />

      {genome.isPending && <LoadingState label="Loading genome" />}
      {genome.error && (
        <ErrorState
          error={genome.error}
          subject={digest}
          apiBase={config.api_base}
          onRetry={() => genome.refetch()}
        />
      )}

      {record && (
        <>
          <GenomeSummaryCard genome={record} aliases={aliasNames} />

          <FhrPanel fhr={record.fhr} />

          <section>
            <h2 className="text-xl font-semibold mb-4">Assets</h2>
            <GenomeAssetTable
              assets={assets.data?.results}
              archives={archives.data?.results}
              loading={assets.isPending}
              error={assets.error}
              onRetry={() => assets.refetch()}
              empty={
                <div className="rg-state rg-state--empty">
                  <p className="mb-4">No assets for this genome yet.</p>
                  <CapabilityGate cap="build">
                    <Link className="rg-btn rg-btn--primary" to={buildHref}>
                      Build the first one
                    </Link>
                  </CapabilityGate>
                </div>
              }
            />
          </section>
        </>
      )}
    </div>
  );
}
