/**
 * The management surface: everything that changes local state but is not a
 * long-running job. Local mode only — the route is registered unconditionally
 * so a bookmark degrades cleanly rather than 404ing, and the capability check
 * decides what it renders.
 */

import { useState } from 'react';
import { Link } from 'react-router-dom';
import { useCapability } from '../hooks/useCapability';
import { SubscriptionPanel } from '../components/manage/SubscriptionPanel';
import { AliasPanel } from '../components/manage/AliasPanel';
import { AssetClassPanel, RecipePanel } from '../components/manage/RecipePanel';
import { InitGenomeModal } from '../components/actions/InitGenomeModal';
import { NotAvailablePage } from './NotAvailablePage';

export function ManagePage() {
  const canSubscribe = useCapability('subscriptions');
  const canWriteAliases = useCapability('aliases_write');
  const canInitGenome = useCapability('genome_init');
  const canBuild = useCapability('build');
  const [initOpen, setInitOpen] = useState(false);

  const anything = canSubscribe || canWriteAliases || canInitGenome || canBuild;
  if (!anything) {
    return (
      <NotAvailablePage
        title="Manage"
        reason="This instance exposes no management operations."
      />
    );
  }

  return (
    <div className="flex flex-col gap-10 max-w-screen-lg mx-auto">
      <header className="flex items-start justify-between gap-4 flex-wrap">
        <h1 className="text-3xl font-bold">Manage</h1>
        <div className="flex gap-2 flex-wrap">
          {canInitGenome && (
            <button
              type="button"
              className="rg-btn rg-btn--primary"
              onClick={() => setInitOpen(true)}
            >
              Initialize genome
            </button>
          )}
          {canBuild && (
            <Link className="rg-btn" to="/build">
              Build asset
            </Link>
          )}
        </div>
      </header>

      <SubscriptionPanel />
      <AliasPanel />
      <RecipePanel />
      <AssetClassPanel />

      <section className="rg-danger-zone">
        <h2 className="text-xl font-semibold">Danger zone</h2>
        <p className="rg-muted text-sm">
          Deleting a genome cascades: its aliases, asset groups and assets all go with it,
          unconditionally. Genome and asset deletion live on their own pages, next to the
          thing being deleted, so the confirmation can name what is about to disappear.
        </p>
        <Link className="rg-link" to="/genomes">
          Go to genomes
        </Link>
      </section>

      <InitGenomeModal isOpen={initOpen} onClose={() => setInitOpen(false)} />
    </div>
  );
}
