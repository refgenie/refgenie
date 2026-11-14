import { useState } from 'react';
import { Copyable } from '../common/Copyable';
import { BridgePullButton } from '../bridge/BridgePullButton';
import { ConnectDialog } from '../bridge/ConnectDialog';
import { useBridge } from '../../hooks/useBridge';
import { formatBytes } from '../../utils/format';
import ga4ghMark from '../../assets/ga4gh-mark.png';

export const SEQCOL_SPEC_URL = 'https://ga4gh.github.io/refget/seqcols/';

export interface GenomeDigestProps {
  digest: string;
  /** The `refgenie pull` command; adds a "Get it with refgenie" row when set. */
  pullCommand?: string | null;
  /** Bytes the pull fetches, shown beside the command so a multi-GB pull is no surprise. */
  pullSize?: number | null;
  /** Whether the genome has a fasta asset a connected local refgenie can pull. */
  hasFastaAsset?: boolean;
}

/**
 * The genome's identity: its GA4GH refget sequence collection digest, in full,
 * with a copy button and a "?" that links to the standard.
 */
export function GenomeDigest({
  digest,
  pullCommand,
  pullSize = null,
  hasFastaAsset = false,
}: GenomeDigestProps) {
  return (
    <div className="rg-genome-digest">
      <div className="rg-genome-digest__row">
        <DigestTerm />
        <Copyable
          value={digest}
          label="Copy digest"
          className="rg-genome-digest__value"
          mono
        >
          {digest}
        </Copyable>
      </div>
      {pullCommand && (
        <GetWithRefgenie
          digest={digest}
          command={pullCommand}
          size={pullSize}
          hasFastaAsset={hasFastaAsset}
        />
      )}
    </div>
  );
}

/** The refgenie row: the command to copy, and a way to do it from the browser. */
function GetWithRefgenie({
  digest,
  command,
  size,
  hasFastaAsset,
}: {
  digest: string;
  command: string;
  size: number | null;
  hasFastaAsset: boolean;
}) {
  const bridge = useBridge();
  const [dialogOpen, setDialogOpen] = useState(false);
  const canConnect =
    bridge.status !== 'connected' &&
    bridge.status !== 'blocked' &&
    bridge.status !== 'self';
  return (
    <div className="rg-genome-digest__row rg-genome-digest__row--get">
      <span className="rg-genome-digest__label">
        <img src="refgenie_bug.svg" alt="" className="rg-genome-digest__mark" />
        <span className="rg-genome-digest__title">Get it with refgenie</span>
      </span>
      <span className="rg-genome-digest__pull">
        <Copyable
          value={command}
          label="Copy refgenie command"
          className="rg-genome-digest__command"
          mono
        >
          {command}
        </Copyable>
        {size !== null && (
          <span className="rg-genome-digest__size" title="Size of the asset once pulled">
            {formatBytes(size)}
          </span>
        )}
      </span>
      {bridge.status === 'connected' && hasFastaAsset && (
        <BridgePullButton genomeDigest={digest} assetGroupName="fasta" />
      )}
      {canConnect && (
        <>
          <button
            type="button"
            className="rg-btn rg-btn--sm rg-genome-digest__connect"
            onClick={() => setDialogOpen(true)}
          >
            Connect local refgenie
          </button>
          <ConnectDialog open={dialogOpen} onClose={() => setDialogOpen(false)} />
        </>
      )}
    </div>
  );
}

/** The ring, "refget digest" and a "?" to the standard; also the Details tab's row label. */
export function DigestTerm({ title = 'refget digest' }: { title?: string }) {
  return (
    <span className="rg-genome-digest__label">
      <img src={ga4ghMark} alt="" className="rg-genome-digest__mark" />
      <span className="rg-genome-digest__title">{title}</span>
      <a
        href={SEQCOL_SPEC_URL}
        className="rg-genome-digest__help"
        target="_blank"
        rel="noreferrer"
        aria-label="What is a refget digest? (opens in a new tab)"
        title="GA4GH refget sequence collections standard"
      >
        ?
      </a>
    </span>
  );
}
