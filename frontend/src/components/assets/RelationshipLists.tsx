import { Link } from 'react-router-dom';
import { DigestChip } from '../common/DigestChip';
import type { AssetResponse } from '../../types/api';

export interface RelationshipListsProps {
  parents: AssetResponse[] | undefined;
  children: AssetResponse[] | undefined;
}

function AssetLinkList({ assets, empty }: { assets: AssetResponse[]; empty: string }) {
  if (assets.length === 0) return <p className="rg-muted text-sm">{empty}</p>;
  return (
    <ul className="flex flex-col gap-2">
      {assets.map((asset) => (
        <li className="flex items-center gap-2 flex-wrap text-sm" key={asset.digest ?? asset.name}>
          {asset.digest ? (
            <Link className="rg-link" to={`/assets/${asset.digest}`}>
              {asset.asset_group_name ? `${asset.asset_group_name}:` : ''}
              {asset.name}
            </Link>
          ) : (
            <span>{asset.name}</span>
          )}
          <DigestChip digest={asset.digest} copyable={false} />
        </li>
      ))}
    </ul>
  );
}

/**
 * Parent and child assets from `/relationships/{digest}?expand=true`. The old
 * dash asset page had none of this, so local mode gains it here.
 */
export function RelationshipLists({ parents, children }: RelationshipListsProps) {
  return (
    <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
      <div>
        <h3 className="text-sm font-semibold mb-2">Parents</h3>
        <AssetLinkList assets={parents ?? []} empty="No parent assets." />
      </div>
      <div>
        <h3 className="text-sm font-semibold mb-2">Children</h3>
        <AssetLinkList assets={children ?? []} empty="No child assets." />
      </div>
    </div>
  );
}
