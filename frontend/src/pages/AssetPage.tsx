import { Link, useNavigate, useParams } from 'react-router-dom';
import { useAsset, useAssetFiles } from '../hooks/queries/useAssets';
import { useRelationshipsExpanded } from '../hooks/queries/useRelationships';
import { useStagedAssets } from '../hooks/queries/useStagedAssets';
import { useUiConfig } from '../hooks/useUiConfig';
import { useApiClient } from '../hooks/useApiClient';
import { useRouteBase } from '../hooks/useRouteBase';
import { Breadcrumbs } from '../components/common/Breadcrumbs';
import { MiniHero } from '../components/layout/MiniHero';
import { ScopedLink } from '../components/common/ScopedLink';
import { BridgePullButton } from '../components/bridge/BridgePullButton';
import { DeleteAssetButton } from '../components/actions/DeleteAssetButton';
import { DescriptionList } from '../components/common/DescriptionList';
import { DigestChip } from '../components/common/DigestChip';
import { ApiLink } from '../components/common/ApiLink';
import { FileSize } from '../components/common/FileSize';
import { Badge } from '../components/common/Badge';
import { ErrorState, LoadingState } from '../components/common/states';
import { ServingModeBadges } from '../components/assets/ServingModeBadges';
import { AssetFileList } from '../components/assets/AssetFileList';
import { RelationshipLists } from '../components/assets/RelationshipLists';
import { StagedAssetTable } from '../components/genomes/StagedAssetTable';

export function AssetPage() {
  const { digest } = useParams<{ digest: string }>();
  const config = useUiConfig();
  const client = useApiClient();
  const routeBase = useRouteBase();
  const navigate = useNavigate();

  const asset = useAsset(digest);
  const files = useAssetFiles(digest);
  const relationships = useRelationshipsExpanded(digest);
  const staged = useStagedAssets({ asset_digest: digest }, { enabled: !!digest });

  // The head renders in every state, so nothing may early-return past it: a
  // loading or missing asset still gets a heading, a breadcrumb out, and (via
  // MiniHero) a correct document title.
  const record = asset.data;
  const label = record
    ? `${record.asset_group_name ? `${record.asset_group_name}:` : ''}${record.name}`
    : 'Asset';
  const builds = (record?.names ?? []).filter((n) => n.build_digest !== null);

  return (
    <div className="flex flex-col gap-8">
      <MiniHero
        title={label}
        // 'Asset' is a placeholder heading while the record loads. Never put it
        // in the tab strip: `false` keeps the bare service name there instead.
        documentTitle={record ? undefined : false}
        breadcrumbs={
          <Breadcrumbs
            items={[
              { label: 'Genomes', to: '/genomes' },
              ...(record?.genome_digest
                ? [{ label: 'Genome', to: `/genomes/${record.genome_digest}` }]
                : []),
              { label },
            ]}
          />
        }
        actions={
          record && (
            <>
              {/* The cross-origin "pull to my refgenie" affordance. It hides
                  itself unless the bridge is connected and the LOCAL instance
                  reports the `pull` capability, so it is inert everywhere
                  else. */}
              {record.genome_digest && record.asset_group_name && (
                <BridgePullButton
                  genomeDigest={record.genome_digest}
                  assetGroupName={record.asset_group_name}
                  assetName={record.name}
                />
              )}
              {record.digest && (
                <DeleteAssetButton
                  digest={record.digest}
                  registryPath={label}
                  seekKeys={record.seek_keys}
                  onDeleted={() =>
                    navigate(
                      record.genome_digest
                        ? `${routeBase}/genomes/${record.genome_digest}`
                        : `${routeBase}/genomes`,
                    )
                  }
                />
              )}
            </>
          )
        }
        // The one detail-page exception to "list pages define the concept":
        // `/assets` was deleted, so nothing above this page defines "asset".
        lede={
          <>
            An asset is one built bundle of files for one genome. Its digest is computed from
            those file contents, so two builds that produce identical output are the same asset
            even when they ran years apart on different machines.
          </>
        }
      />

      {asset.isPending && <LoadingState label="Loading asset" />}
      {asset.error && (
        <ErrorState error={asset.error} subject={digest} apiBase={config.api_base} />
      )}

      {record && (
        <>
          <DescriptionList
            items={[
              { term: 'Digest', value: <DigestChip digest={record.digest} length={32} /> },
              { term: 'Description', value: record.description ?? '' },
              { term: 'Size', value: <FileSize bytes={record.size} /> },
              {
                term: 'Names',
                value: record.names?.length ? record.names.map((n) => n.name).join(', ') : '',
              },
              {
                term: 'Serving modes',
                value: <ServingModeBadges modes={record.serving_modes} />,
              },
              { term: 'Asset class', value: record.asset_class_name ?? '' },
              {
                term: 'Asset group',
                value:
                  record.asset_group_id !== null && record.asset_group_id !== undefined ? (
                    <ScopedLink
                      className="rg-link"
                      to={`/asset-groups/${record.asset_group_id}`}
                    >
                      {record.asset_group_name ?? `#${record.asset_group_id}`}
                    </ScopedLink>
                  ) : (
                    ''
                  ),
              },
              {
                term: 'Recipe',
                // Guarded: the Jinja page linked /recipes/None when this was null.
                value:
                  record.recipe_id !== null && record.recipe_id !== undefined ? (
                    <Link className="rg-link" to={`/recipes/${record.recipe_id}`}>
                      #{record.recipe_id}
                    </Link>
                  ) : (
                    ''
                  ),
              },
            ]}
          />

          <section>
            <h2 className="text-xl font-semibold mb-4">Seek keys</h2>
            {record.seek_keys && record.seek_keys.length > 0 ? (
              <ul className="flex flex-col gap-2">
                {record.seek_keys.map((key) => (
                  <li className="flex flex-wrap items-center gap-2 text-sm" key={key.name}>
                    <strong>{key.name}</strong>
                    <Badge>{key.type}</Badge>
                    <code className="rg-code rg-code--inline">{key.value}</code>
                    {key.size !== null && key.size !== undefined && (
                      <span className="rg-muted">
                        (<FileSize bytes={key.size} />)
                      </span>
                    )}
                    {key.description && <span className="rg-muted">— {key.description}</span>}
                  </li>
                ))}
              </ul>
            ) : (
              <p className="rg-muted text-sm">This asset declares no seek keys.</p>
            )}
          </section>

          {/* Provenance is per-name, not per-asset: several builds can produce this
              one content digest, and each records itself against the name it was
              built under, so it gets its own section rather than generic seek
              keys above. */}
          <section>
            <h2 className="text-xl font-semibold mb-4">Builds</h2>
            {builds.length > 0 ? (
              <ul className="flex flex-col gap-4">
                {builds.map((build) => (
                  <li className="flex flex-col gap-1 text-sm" key={build.name}>
                    <div className="flex flex-wrap items-center gap-2">
                      <strong>{build.name}</strong>
                      {build.is_default && <Badge>default</Badge>}
                      <DigestChip digest={build.build_digest} length={16} />
                    </div>
                    <div className="rg-muted flex flex-wrap gap-4">
                      {build.build_timestamp && <span>{build.build_timestamp}</span>}
                      {build.refgenie_version && <span>refgenie {build.refgenie_version}</span>}
                      {build.docker_image && <span>{build.docker_image}</span>}
                    </div>
                    {build.docker_image_digest && (
                      <code className="rg-code rg-code--inline">{build.docker_image_digest}</code>
                    )}
                    {build.inputs && Object.keys(build.inputs).length > 0 && (
                      <pre className="rg-code text-xs overflow-x-auto">
                        {JSON.stringify(build.inputs, null, 2)}
                      </pre>
                    )}
                  </li>
                ))}
              </ul>
            ) : (
              <p className="rg-muted text-sm">
                No build is recorded for this asset's names. Assets that were pulled, or
                built before provenance moved onto the name, show nothing here.
              </p>
            )}
          </section>

          <section>
            <h2 className="text-xl font-semibold mb-4">Files</h2>
            {files.isPending ? (
              <LoadingState label="Loading files" />
            ) : files.error ? (
              <p className="rg-muted text-sm">
                This asset is not staged for file-level serving.
              </p>
            ) : (
              <AssetFileList assetDigest={record.digest ?? ''} files={files.data?.files} />
            )}
          </section>

          <section>
            <h2 className="text-xl font-semibold mb-4">Relationships</h2>
            {relationships.isPending ? (
              <LoadingState label="Loading relationships" />
            ) : relationships.error ? (
              <ErrorState error={relationships.error} apiBase={config.api_base} />
            ) : (
              <RelationshipLists
                parents={relationships.data?.parents}
                children={relationships.data?.children}
              />
            )}
          </section>

          <section>
            <h2 className="text-xl font-semibold mb-4">Staging</h2>
            <StagedAssetTable
              staged={staged.data?.items}
              loading={staged.isPending}
              error={staged.error}
            />
          </section>

          <section>
            <h2 className="text-xl font-semibold mb-4">API endpoints</h2>
            <ul className="flex flex-col gap-1 text-sm">
              <li>
                <ApiLink href={client.url(`/assets/${record.digest}`)}>Asset JSON</ApiLink>
              </li>
              <li>
                <ApiLink href={client.url(`/assets/${record.digest}/files`)}>File list JSON</ApiLink>
              </li>
              <li>
                <ApiLink href={client.url(`/relationships/${record.digest}`, { expand: true })}>
                  Relationships JSON
                </ApiLink>
              </li>
            </ul>
          </section>
        </>
      )}
    </div>
  );
}
