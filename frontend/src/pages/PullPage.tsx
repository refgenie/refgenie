/**
 * `/pull?server=&genome=&asset_group=`: a remote UI's hand-off to the local
 * dash. Shows a prefilled pull for the user to confirm; it must never run from
 * the URL alone, which would be a CSRF hole.
 */

import { Link, Navigate, useNavigate, useSearchParams } from 'react-router-dom';
import { useCapability } from '../hooks/useCapability';
import { usePullAction } from '../hooks/usePullAction';
import { MiniHero } from '../components/layout/MiniHero';
import { DescriptionList } from '../components/common/DescriptionList';
import { NotAvailablePage } from './NotAvailablePage';

export function PullPage() {
  const canPull = useCapability('pull');
  const [params] = useSearchParams();
  const navigate = useNavigate();
  const { pull, isPending } = usePullAction();

  // `genome` is a digest.
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

  if (!genome || !assetGroup) {
    return <Navigate to="/remote" replace />;
  }

  const submit = async () => {
    const submitted = await pull({
      genomeDigest: genome,
      assetGroupName: assetGroup,
      assetName: asset,
      serverUrl: server,
      force: false,
    });
    if (submitted) navigate('/jobs');
  };

  return (
    <div className="flex flex-col gap-8">
      <MiniHero
        title="Pull asset"
        documentTitle="Pull"
        lede={
          <>
            A remote page handed off this pull request. Nothing has been submitted yet. Review
            it, then pull.
          </>
        }
      />

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
