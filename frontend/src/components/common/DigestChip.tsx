import { Copyable } from './Copyable';
import { formatDigest } from '../../utils/format';

export interface DigestChipProps {
  digest: string | null | undefined;
  /** Hide the copy button in dense table cells. */
  copyable?: boolean;
  length?: number;
  /** The copy button's accessible name. */
  label?: string;
  /** `plate`: the white bordered look of the header's refget digest, for detail lists. */
  variant?: 'chip' | 'plate';
}

/** Digests always render through this: monospace, truncated, full value in `title`. */
export function DigestChip({
  digest,
  copyable = true,
  length = 12,
  label = 'Copy digest',
  variant = 'chip',
}: DigestChipProps) {
  if (!digest) return <span className="rg-muted">NA</span>;
  const text = (
    <code className="rg-digest__text" title={digest}>
      {formatDigest(digest, length)}
    </code>
  );
  if (!copyable) return <span className="rg-digest">{text}</span>;
  return (
    <Copyable
      value={digest}
      label={label}
      className={variant === 'plate' ? 'rg-digest rg-digest--plate' : 'rg-digest'}
    >
      {text}
    </Copyable>
  );
}
