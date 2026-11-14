import { Copyable } from '../common/Copyable';
import { ExternalLink } from '../common/ExternalLink';
import { assetRequestUrl } from '../../services/projectLinks';

export interface NoAssetsPanelProps {
  digest: string;
  alias?: string;
  /** A short real range of this genome's first sequence, e.g. `chr1:0-100`. */
  region?: string;
}

function Command({ value }: { value: string }) {
  return (
    <Copyable value={value} label={`Copy: ${value}`} className="rg-command">
      <code className="rg-code rg-code--inline">{value}</code>
    </Copyable>
  );
}

/** What a visitor can do on a genome that has no built assets yet. */
export function NoAssetsPanel({ digest, alias, region }: NoAssetsPanelProps) {
  return (
    <div className="rg-state rg-state--empty">
      <h2>No assets yet</h2>
      <p className="mb-4">
        Nobody has built assets for this genome on this server yet. Assets are files
        made from the genome, such as aligner indexes and annotations.
      </p>
      {alias && region && (
        <div className="mb-4">
          <p className="mb-2">Only need the sequence?</p>
          <Command value={`refgenie pull --init -g ${alias}`} />
          <Command value={`refgenie getseq -g ${alias} -l ${region}`} />
          <p className="mt-2 text-sm">
            <code>--init</code> registers the genome and downloads no files.
          </p>
        </div>
      )}
      <p>
        Need an asset built?{' '}
        <ExternalLink href={assetRequestUrl(alias, digest)}>Request assets</ExternalLink>{' '}
        (a GitHub issue on refgenie-registry, filled in with this genome).
      </p>
    </div>
  );
}
