import { useMemo } from 'react';
import { Link, useParams, useSearchParams } from 'react-router-dom';
import { useGenome } from '../hooks/queries/useGenomes';
import { useAssets } from '../hooks/queries/useAssets';
import { useAliases } from '../hooks/queries/useAliases';
import { useArchives } from '../hooks/queries/useServerInfo';
import { useGenomeSequences } from '../hooks/queries/useGenomeSequences';
import { useCapability } from '../hooks/useCapability';
import { useFastaDownload } from '../hooks/useFastaDownload';
import { useCompareCandidates } from '../hooks/useCompareCandidates';
import { useUiConfig } from '../hooks/useUiConfig';
import { preferredAlias } from '../utils/genomes';
import { namingStyle } from '../utils/namingStyle';
import { formatDigest } from '../utils/format';
import { fastaDownloadAvailable, fastaPullSize, isGenomeTab } from '../utils/genomePage';
import type { GenomeTab } from '../utils/genomePage';
import { ApiError } from '../services/http';
import { fastaDownloadEnabled } from '../services/config';
import { Breadcrumbs } from '../components/common/Breadcrumbs';
import { MiniHero } from '../components/layout/MiniHero';
import { CapabilityGate } from '../components/common/CapabilityGate';
import { DeleteGenomeButton } from '../components/actions/DeleteGenomeButton';
import { ErrorState, LoadingState } from '../components/common/states';
import { TabPanel, Tabs } from '../components/common/Tabs';
import type { TabItem } from '../components/common/Tabs';
import { GenomeDigest } from '../components/genomes/GenomeDigest';
import { GenomeHeader } from '../components/genomes/GenomeHeader';
import { refgeniePullCommand } from '../utils/codeSnippets';
import { GenomeSummaryCard } from '../components/genomes/GenomeSummaryCard';
import { GenomeAliases } from '../components/genomes/GenomeAliases';
import { FhrPanel } from '../components/genomes/FhrPanel';
import { GenomeAssetTable } from '../components/genomes/GenomeAssetTable';
import { GenomeDownloadPanel } from '../components/genomes/GenomeDownload';
import {
  ComparePanel,
  ExtractPanel,
  SequencesPanel,
} from '../components/sequences/GenomeSequencePanels';
import { CollectionDigests } from '../components/sequences/CollectionDigests';
import { isExtractable } from '../services/regionExtract';
import { assetRequestUrl } from '../services/projectLinks';
import { ExternalLink } from '../components/common/ExternalLink';
import { NoAssetsPanel } from '../components/genomes/NoAssetsPanel';

const NO_SEQUENCES = 'No sequence data for this genome';
const NO_SAME_SPECIES = 'No other genomes of this species';
const NO_STORE = 'This genome has no refget store';
const NO_READABLE = 'No sequence here can be read by range';

export function GenomePage() {
  const { digest } = useParams<{ digest: string }>();
  const config = useUiConfig();
  const canReadArchives = useCapability('archives');
  const canBuild = useCapability('build');

  const genome = useGenome(digest);
  const aliases = useAliases({ genome_digest: digest, pageSize: 100 });
  const assets = useAssets(
    { genome_digest: digest, pageSize: 200 },
    { enabled: !!digest },
  );
  const archives = useArchives(
    { genome_digest: digest },
    { enabled: canReadArchives && !!digest },
  );
  const sequences = useGenomeSequences(digest);
  const [params, setParams] = useSearchParams();

  // No early return: the head renders in every state.
  const record = genome.data;
  const aliasNames = aliases.data?.results.map((alias) => alias.name) ?? [];
  // A new genome may have no alias; fall back to a shortened digest.
  const primaryAlias = preferredAlias(aliasNames);
  const title = record ? (primaryAlias ?? formatDigest(record.digest)) : 'Genome';

  const storeUrl = record?.refget_store_url ?? null;
  // The browser reads the store itself only when the build turns that on.
  const canDownload = !!record && fastaDownloadAvailable(record);
  const fasta = useFastaDownload({
    storeUrl: canDownload ? record.refget_store_url : null,
    digest: digest ?? '',
    fileStem: primaryAlias ?? digest ?? 'genome',
  });

  // A 404 means the genome has no sequence data here: a normal answer.
  const noSequences = sequences.error instanceof ApiError && sequences.error.isNotFound;
  const hasSequences = sequences.data
    ? true
    : noSequences
      ? false
      : sequences.error
        ? true
        : null;
  const sequenceStats = sequences.data
    ? {
        count: sequences.data.sequences.length,
        bases: sequences.data.sequences.reduce((sum, s) => sum + s.length, 0),
      }
    : hasSequences === false || sequences.error
      ? (false as const)
      : null;
  const naming = useMemo(
    () => (sequences.data ? namingStyle(sequences.data.sequences) : null),
    [sequences.data],
  );
  const assetCount = assets.data ? assets.data.results.length : null;
  const hasFastaAsset = !!assets.data?.results.some(
    (asset) => asset.asset_class_name === 'fasta',
  );
  // A local dash already has the genome; the command is for visitors.
  const pullCommand =
    config.mode === 'local' || !assets.data
      ? null
      : refgeniePullCommand(primaryAlias, hasFastaAsset);
  const candidates = useCompareCandidates(record);
  const noCandidates = candidates !== null && candidates.length === 0;
  // The in-browser FASTA comparison keeps Compare open without other genomes.
  const compareFile = fastaDownloadEnabled();
  const compareOff =
    hasSequences === false
      ? NO_SEQUENCES
      : noCandidates && !compareFile
        ? NO_SAME_SPECIES
        : undefined;
  const readable = sequences.data
    ? sequences.data.sequences.some((row) => isExtractable(row.alphabet))
    : null;
  const downloadOff = record && !storeUrl ? NO_STORE : undefined;
  const extractOff =
    downloadOff ??
    (hasSequences === false ? NO_SEQUENCES : readable === false ? NO_READABLE : undefined);

  const tabs: TabItem<GenomeTab>[] = [
    {
      id: 'sequences',
      label: 'Sequences',
      count:
        hasSequences === false ? undefined : sequenceStats ? sequenceStats.count : null,
      disabled: hasSequences === false,
      title: hasSequences === false ? NO_SEQUENCES : undefined,
    },
    { id: 'extract', label: 'Extract', disabled: !!extractOff, title: extractOff },
    { id: 'download', label: 'Download', disabled: !!downloadOff, title: downloadOff },
    { id: 'assets', label: 'Assets', count: assetCount },
    {
      id: 'compare',
      label: 'Compare',
      count: compareOff || noCandidates ? undefined : candidates ? candidates.length : null,
      disabled: !!compareOff,
      title: compareOff,
    },
    { id: 'details', label: 'Details' },
  ];

  const requested = params.get('tab');
  const requestedTab =
    isGenomeTab(requested) && !tabs.find((tab) => tab.id === requested)?.disabled
      ? requested
      : null;
  const activeTab = requestedTab ?? 'assets';
  const selectTab = (id: GenomeTab) => {
    const next = new URLSearchParams(params);
    next.set('tab', id);
    setParams(next);
  };

  // A short real range for the "Only need the sequence?" command.
  const firstSequence = sequences.data?.sequences[0];
  const region = firstSequence
    ? `${firstSequence.name}:0-${Math.min(100, firstSequence.length)}`
    : undefined;

  // The build form needs an alias; the digest is only for the job target.
  const buildHref =
    `/build?genome=${encodeURIComponent(aliasNames[0] ?? '')}` +
    `&genome_digest=${encodeURIComponent(record?.digest ?? '')}`;

  return (
    <div className="flex flex-col gap-6">
      <MiniHero
        title={title}
        // `false` while loading keeps the placeholder 'Genome' out of the tab title.
        documentTitle={record ? undefined : false}
        breadcrumbs={
          <Breadcrumbs
            items={[{ label: 'Genomes', to: '/genomes' }, { label: title }]}
          />
        }
        actions={
          record && (
            <>
              <CapabilityGate cap="build">
                <Link className="rg-btn" to={buildHref}>
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
        aside={
          record ? (
            <GenomeDigest
              digest={record.digest}
              pullCommand={pullCommand}
              pullSize={hasFastaAsset ? fastaPullSize(assets.data?.results) : null}
              hasFastaAsset={hasFastaAsset}
            />
          ) : undefined
        }
      >
        {!genome.error && (
          <GenomeHeader genome={record} naming={naming} />
        )}
      </MiniHero>

      {genome.error && (
        <ErrorState
          error={genome.error}
          subject={digest}
          apiBase={config.api_base}
          onRetry={() => genome.refetch()}
        />
      )}

      {!genome.error && (
        <div>
          <Tabs
            idPrefix="genome"
            label="Genome sections"
            tabs={tabs}
            selected={record ? activeTab : null}
            onSelect={selectTab}
          />
          {record ? (
            <TabPanel idPrefix="genome" id={activeTab}>
              {activeTab === 'assets' && (
                <GenomeAssetTable
                  assets={assets.data?.results}
                  archives={archives.data?.results}
                  loading={assets.isPending}
                  error={assets.error}
                  onRetry={() => assets.refetch()}
                  empty={
                    canBuild ? (
                      <div className="rg-state rg-state--empty">
                        <p className="mb-4">No assets for this genome yet.</p>
                        <Link className="rg-btn rg-btn--primary" to={buildHref}>
                          Build the first one
                        </Link>
                      </div>
                    ) : (
                      <NoAssetsPanel
                        digest={record.digest}
                        alias={primaryAlias}
                        region={region}
                      />
                    )
                  }
                />
              )}
              {activeTab === 'assets' && !canBuild && assetCount !== null && assetCount > 0 && (
                <p className="mt-4 text-sm">
                  Need an asset that is not listed?{' '}
                  <ExternalLink href={assetRequestUrl(primaryAlias, record.digest)}>
                    Request it
                  </ExternalLink>
                </p>
              )}
              {activeTab === 'sequences' && (
                <SequencesPanel digest={record.digest} canExtract={!extractOff} />
              )}
              {activeTab === 'extract' && storeUrl && (
                <ExtractPanel
                  digest={record.digest}
                  storeUrl={storeUrl}
                  browserStoreUrl={canDownload ? storeUrl : null}
                />
              )}
              {activeTab === 'download' && storeUrl && (
                <GenomeDownloadPanel
                  digest={record.digest}
                  storeUrl={storeUrl}
                  alias={primaryAlias}
                  download={canDownload ? fasta : null}
                  hasFastaAsset={hasFastaAsset}
                />
              )}
              {activeTab === 'compare' && (
                <ComparePanel
                  digest={record.digest}
                  name={primaryAlias}
                  candidates={candidates}
                  canCompareFile={compareFile}
                />
              )}
              {activeTab === 'details' && (
                <>
                  <section className="rg-flat-section">
                    <h2>Genome</h2>
                    <GenomeSummaryCard genome={record} />
                  </section>
                  {sequences.data && (
                    <CollectionDigests
                      digest={record.digest}
                      level1={sequences.data.level1}
                      storeUrl={record.refget_store_url}
                    />
                  )}
                  <GenomeAliases
                    aliases={aliasNames}
                    collectionAliases={sequences.data?.collection_aliases}
                  />
                  <FhrPanel fhr={record.fhr} flat />
                </>
              )}
            </TabPanel>
          ) : (
            <div className="rg-tabs__panel">
              <LoadingState rows={8} label="Loading" />
            </div>
          )}
        </div>
      )}
    </div>
  );
}
