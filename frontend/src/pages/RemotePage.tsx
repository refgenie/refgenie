import { useMemo, useState } from 'react';
import { useCapability } from '../hooks/useCapability';
import {
  useRemoteAssets,
  useRemoteGenomes,
  useRemoteServers,
} from '../hooks/queries/useRemote';
import { useGenomes } from '../hooks/queries/useGenomes';
import { useAssets } from '../hooks/queries/useAssets';
import { RemoteServerPanel } from '../components/remote/RemoteServerPanel';
import { RemoteGenomeTable } from '../components/remote/RemoteGenomeTable';
import { RemoteAssetTable } from '../components/remote/RemoteAssetTable';
import { NotAvailablePage } from './NotAvailablePage';

export function RemotePage() {
  const enabled = useCapability('remote_browse');
  const [selectedServer, setSelectedServer] = useState<string | undefined>(undefined);
  const [expanded, setExpanded] = useState<string | undefined>(undefined);

  const servers = useRemoteServers({ enabled });
  const remoteGenomes = useRemoteGenomes(selectedServer, { enabled });
  const remoteAssets = useRemoteAssets(selectedServer, expanded, { enabled: enabled && !!expanded });

  // "Already local?" is a cross-reference against the local catalog.
  const localGenomes = useGenomes({ limit: 1000 });
  const localAssets = useAssets({ limit: 1000 }, { enabled });

  const subscribedCount = (servers.data?.servers ?? []).filter(
    (server) => server.subscribed,
  ).length;

  const localGenomeDigests = useMemo(
    () => new Set((localGenomes.data?.items ?? []).map((genome) => genome.digest)),
    [localGenomes.data],
  );
  const localAssetDigests = useMemo(
    () =>
      new Set(
        (localAssets.data?.items ?? [])
          .map((asset) => asset.digest)
          .filter((digest): digest is string => !!digest),
      ),
    [localAssets.data],
  );

  // The route is registered unconditionally so a bookmarked URL degrades
  // cleanly rather than 404ing.
  if (!enabled) {
    return (
      <NotAvailablePage
        title="Remote browse"
        reason="This instance does not expose remote browsing (capability `remote_browse` is off)."
      />
    );
  }

  return (
    <div className="flex flex-col gap-8">
      <h1 className="text-3xl font-bold">Remote assets</h1>

      <section>
        <h2 className="text-xl font-semibold mb-4">Servers</h2>
        <RemoteServerPanel
          servers={servers.data?.servers}
          selected={selectedServer}
          onSelect={(url) => {
            setSelectedServer(url);
            setExpanded(undefined);
          }}
        />
      </section>

      <section>
        <h2 className="text-xl font-semibold mb-4">Genomes</h2>
        <RemoteGenomeTable
          genomes={remoteGenomes.data}
          localDigests={localGenomeDigests}
          expanded={expanded}
          onToggle={(digest) => setExpanded((current) => (current === digest ? undefined : digest))}
          loading={remoteGenomes.isPending}
          error={remoteGenomes.error}
        />
      </section>

      {expanded && (
        <section>
          <h2 className="text-xl font-semibold mb-4">Assets</h2>
          <RemoteAssetTable
            assets={remoteAssets.data}
            localDigests={localAssetDigests}
            loading={remoteAssets.isPending}
            error={remoteAssets.error}
            // Only meaningful when there is more than one server to move to.
            onPickServer={
              subscribedCount > 1 ? () => setSelectedServer(undefined) : undefined
            }
          />
        </section>
      )}
    </div>
  );
}
