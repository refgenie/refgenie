/**
 * `/pull` — the localhost-bridge deep-link hand-off.
 *
 * The public SPA (and any other remote UI) hands the browser off to
 * `http://localhost:{port}/pull?server={remoteApiBase}&genome={digest}&asset_group={name}`
 * (documented next to `SPA_CLIENT_ROUTES` in `refgenie/server/const.py`).
 * This page parses those params and renders a prefilled, UNSUBMITTED pull
 * confirmation. It must never auto-execute from the URL: that would recreate
 * the CSRF hole the X-Refgenie-Action header closes.
 *
 * Incomplete params degrade to the remote-browse view; a missing `pull`
 * capability degrades to the standard not-available page.
 */

import { Link, Navigate, useNavigate, useSearchParams } from 'react-router-dom';
import { useCapability } from '../hooks/useCapability';
import { usePullAction } from '../hooks/usePullAction';
import { DescriptionList } from '../components/common/DescriptionList';
import { NotAvailablePage } from './NotAvailablePage';

export function PullPage() {
  const canPull = useCapability('pull');
  const [params] = useSearchParams();
  const navigate = useNavigate();
  const { pull, isPending } = usePullAction();

  // The documented hand-off params. `genome` is a digest in the contract.
  const server = params.get('server');
  const genome = params.get('genome');
  const assetGroup = params.get('asset_group');
  const asset = params.get('asset');

  if (!canPull) {
    return (
      <NotAvailablePage
        title="Pull"
        reason="This instance cannot pull assets (capability `pull` is off)."
      />
    );
  }

  // Without a target there is nothing to confirm; browsing is the next best.
  if (!genome || !assetGroup) {
    return <Navigate to="/remote" replace />;
  }

  const submit = async () => {
    const ref = await pull({
      genomeDigest: genome,
      assetGroupName: assetGroup,
      assetName: asset,
      serverUrl: server,
      force: false,
    });
    // The console is where the queued job's progress lives.
    if (ref) navigate('/jobs');
  };

  return (
    <div className="flex flex-col gap-8">
      <h1 className="text-3xl font-bold">Pull asset</h1>
      <p className="rg-muted">
        A remote page handed off this pull request. Nothing has been submitted
        yet — review it, then pull.
      </p>

      <DescriptionList
        hideEmpty
        items={[
          { term: 'Genome', value: <code className="rg-code rg-code--inline">{genome}</code> },
          { term: 'Asset group', value: assetGroup },
          { term: 'Asset', value: asset },
          { term: 'Server', value: server },
        ]}
      />

      <div className="flex items-center gap-4">
        <button
          type="button"
          className="rg-btn rg-btn--primary"
          disabled={isPending}
          onClick={() => void submit()}
        >
          Pull
        </button>
        <Link className="rg-link" to="/remote">
          Browse remote assets instead
        </Link>
      </div>
    </div>
  );
}
