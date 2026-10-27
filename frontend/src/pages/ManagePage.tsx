/**
 * The management surface: everything that changes local state but is not a
 * long-running job. Local mode only — the route is registered unconditionally
 * so a bookmark degrades cleanly rather than 404ing, and the capability check
 * decides what it renders.
 */

import { useState } from 'react';
import { Link } from 'react-router-dom';
import { useCapability } from '../hooks/useCapability';
import { MiniHero } from '../components/layout/MiniHero';
import { SubscriptionPanel } from '../components/manage/SubscriptionPanel';
import { DataChannelPanel } from '../components/manage/DataChannelPanel';
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
    // Same rule as BuildPage: the measure cap goes on the panels, not around
    // the MiniHero, whose full-bleed margin is computed from its parent's
    // width. This page only escaped that bug because `mx-auto` happened to
    // centre the capped wrapper, which is the one arrangement where half the
    // parent still lands where the head expects.
    <div className="flex flex-col gap-8">
      <MiniHero
        title="Manage"
        lede={
          <>
            Settings for this refgenie install: which servers it pulls from, which names point
            at which genomes, and which recipes and asset classes it has registered. Registering
            a new recipe or asset class is still a command-line operation.
          </>
        }
        actions={
          <>
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
          </>
        }
      />

      <div className="flex flex-col gap-8 max-w-screen-lg mx-auto w-full">
        <SubscriptionPanel />
        <DataChannelPanel />
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
      </div>

      <InitGenomeModal isOpen={initOpen} onClose={() => setInitOpen(false)} />
    </div>
  );
}
