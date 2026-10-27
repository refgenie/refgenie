import { CopyButton } from './CopyButton';
import { formatDigest } from '../../utils/format';

export interface DigestChipProps {
  digest: string | null | undefined;
  /** Hide the copy button in dense table cells. */
  copyable?: boolean;
  length?: number;
}

/** Digests always render through this: monospace, truncated, full value in `title`. */
export function DigestChip({ digest, copyable = true, length = 12 }: DigestChipProps) {
  if (!digest) return <span className="rg-muted">NA</span>;
  return (
    <span className="rg-digest">
      <code className="rg-digest__text" title={digest}>
        {formatDigest(digest, length)}
      </code>
      {copyable && <CopyButton value={digest} label="Copy digest" />}
    </span>
  );
}
