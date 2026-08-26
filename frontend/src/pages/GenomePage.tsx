import { Link, useParams } from 'react-router-dom';
import { useGenome } from '../hooks/queries/useGenomes';
import { useAssets } from '../hooks/queries/useAssets';
import { useAliases } from '../hooks/queries/useAliases';
import { useArchives } from '../hooks/queries/useServerInfo';
import { useCapability } from '../hooks/useCapability';
import { useUiConfig } from '../hooks/useUiConfig';
import { Breadcrumbs } from '../components/common/Breadcrumbs';
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
  const aliases = useAliases({ genome_digest: digest, limit: 100 });
  const assets = useAssets({ genome_digest: digest, limit: 200 }, { enabled: !!digest });
  const archives = useArchives({ genome_digest: digest }, { enabled: canReadArchives && !!digest });

  if (genome.isPending) return <LoadingState label="Loading genome" />;
  if (genome.error) {
    return <ErrorState error={genome.error} subject={digest} apiBase={config.api_base} />;
  }
  if (!genome.data) return null;

  const aliasNames = aliases.data?.items.map((alias) => alias.name) ?? [];
  const title = aliasNames[0] ?? genome.data.digest;
  // The build form wants the ALIAS: preflight resolves it with `alias.resolve`,
  // which does not accept a digest. The digest rides along for the job target.
  const buildHref =
    `/build?genome=${encodeURIComponent(aliasNames[0] ?? '')}` +
    `&genome_digest=${encodeURIComponent(genome.data.digest)}`;

  return (
    <div className="flex flex-col gap-8">
      <Breadcrumbs items={[{ label: 'Genomes', to: '/genomes' }, { label: title }]} />

      <header className="flex items-start justify-between gap-4 flex-wrap">
        <h1 className="text-3xl font-bold">{title}</h1>
        <div className="flex gap-2 flex-wrap">
          <CapabilityGate cap="build">
            <Link
              className="rg-btn rg-btn--primary"
              to={buildHref}
            >
              Build asset
            </Link>
          </CapabilityGate>
          <DeleteGenomeButton
            digest={genome.data.digest}
            primaryAlias={aliasNames[0] ?? genome.data.digest}
            aliasCount={aliasNames.length}
            assetCount={assets.data?.items.length ?? 0}
          />
        </div>
      </header>

      <GenomeSummaryCard genome={genome.data} aliases={aliasNames} />

      <FhrPanel fhr={genome.data.fhr} />

      <section>
        <h2 className="text-xl font-semibold mb-4">Assets</h2>
        <GenomeAssetTable
          assets={assets.data?.items}
          archives={archives.data?.items}
          loading={assets.isPending}
          error={assets.error}
          onRetry={() => assets.refetch()}
          empty={
            <div className="rg-state rg-state--empty">
              <p className="mb-4">No assets for this genome yet.</p>
              <CapabilityGate cap="build">
                <Link
                  className="rg-btn rg-btn--primary"
                  to={buildHref}
                >
                  Build the first one
                </Link>
              </CapabilityGate>
            </div>
          }
        />
      </section>
    </div>
  );
}
