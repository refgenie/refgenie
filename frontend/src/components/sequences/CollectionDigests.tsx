import { useState } from 'react';
import { Link } from 'react-router-dom';
import { BaseModal } from '../common/BaseModal';
import { DigestChip } from '../common/DigestChip';
import { ExternalLink } from '../common/ExternalLink';
import { DigestTerm } from '../genomes/GenomeDigest';
import { useAliases } from '../../hooks/queries/useAliases';
import { useApiClient } from '../../hooks/useApiClient';
import { useResource } from '../../hooks/useResource';
import { useUiConfig } from '../../hooks/useUiConfig';
import {
  attributeUrl,
  collectionUrl,
  findCollectionsSharing,
  seqcolBase,
  serviceForStore,
} from '../../services/resources/seqcol';
import { formatDigest } from '../../utils/format';
import { preferredAlias } from '../../utils/genomes';

export interface CollectionDigestsProps {
  /** The genome's top-level sequence collection digest. */
  digest: string;
  /** Level-1 attribute digests: names, lengths, sequences, ... */
  level1: Record<string, string>;
  /** The genome's RefgetStore, to pick the seqcol API that serves it. */
  storeUrl?: string | null;
}

/**
 * The Details tab's refget section: the top-level digest, then each level-1
 * attribute digest, each with a link to the seqcol API and (for attributes) a
 * search for other genomes that share it.
 */
export function CollectionDigests({
  digest,
  level1,
  storeUrl,
}: CollectionDigestsProps) {
  const client = useApiClient();
  const services = useUiConfig().seqcol_services ?? [];
  const service = serviceForStore(services, storeUrl);
  const base = service ? seqcolBase(client.baseUrl, service) : null;

  return (
    <section className="rg-flat-section rg-refget">
      <h2 className="rg-refget__head">
        <DigestTerm title="refget sequence collection" />
      </h2>
      <table className="rg-refget__table">
        <tbody>
          <tr>
            <th scope="row" className="rg-refget__term">
              Top-level collection digest
            </th>
            <td>
              <DigestChip digest={digest} length={32} variant="plate" />
            </td>
            <td className="rg-refget__links">
              {base && <ApiLink href={collectionUrl(base, digest)} />}
            </td>
          </tr>
          <tr>
            <th colSpan={3} scope="rowgroup" className="rg-refget__group">
              Level-1 attribute digests
            </th>
          </tr>
          {Object.entries(level1).map(([attribute, value]) => (
            <AttributeRow
              key={attribute}
              attribute={attribute}
              value={value}
              self={digest}
              base={base}
            />
          ))}
        </tbody>
      </table>
    </section>
  );
}

function ApiLink({ href }: { href: string }) {
  return (
    <ExternalLink href={href} icon="code" hint="opens the seqcol API JSON">
      API
    </ExternalLink>
  );
}

function AttributeRow({
  attribute,
  value,
  self,
  base,
}: {
  attribute: string;
  value: string;
  self: string;
  base: string | null;
}) {
  const [open, setOpen] = useState(false);
  return (
    <tr className="rg-refget__attr">
      <th scope="row" className="rg-refget__term rg-refget__term--attr">
        <code>{attribute}</code>
      </th>
      <td>
        <DigestChip
          digest={value}
          length={32}
          label={`Copy ${attribute} digest`}
          variant="plate"
        />
      </td>
      <td className="rg-refget__links">
        {base && <ApiLink href={attributeUrl(base, attribute, value)} />}
        <button
          type="button"
          className="rg-btn rg-btn--sm rg-btn--bare"
          aria-haspopup="dialog"
          onClick={() => setOpen(true)}
        >
          Genomes that share this
        </button>
        <BaseModal
          isOpen={open}
          onClose={() => setOpen(false)}
          title={`Genomes with the same ${attribute}`}
        >
          <div className="text-sm">
            <SharedGenomes attribute={attribute} value={value} self={self} />
          </div>
        </BaseModal>
      </td>
    </tr>
  );
}

const MAX_SHOWN = 25;

function SharedGenomes({
  attribute,
  value,
  self,
}: {
  attribute: string;
  value: string;
  self: string;
}) {
  const client = useApiClient();
  const services = useUiConfig().seqcol_services ?? [];
  const found = useResource(['seqcol-share', attribute, value], ({ signal }) =>
    findCollectionsSharing(client.baseUrl, services, attribute, value, { signal }),
  );
  if (found.error)
    return <p className="rg-muted">Could not search the sequence collections.</p>;
  if (!found.data) return <p className="rg-muted">Searching…</p>;
  const others = found.data.filter((d) => d !== self);
  if (others.length === 0) {
    return (
      <p className="rg-muted">
        No other genome here has the same <code>{attribute}</code>.
      </p>
    );
  }
  return (
    <div>
      <p className="rg-muted">
        {others.length.toLocaleString()} other{' '}
        {others.length === 1 ? 'genome has' : 'genomes have'} the same{' '}
        <code>{attribute}</code>:
      </p>
      <ul className="rg-refget__shared-list">
        {others.slice(0, MAX_SHOWN).map((d) => (
          <li key={d}>
            <SharedGenomeLink digest={d} />
          </li>
        ))}
      </ul>
      {others.length > MAX_SHOWN && (
        <p className="rg-muted">
          and {(others.length - MAX_SHOWN).toLocaleString()} more.
        </p>
      )}
    </div>
  );
}

function SharedGenomeLink({ digest }: { digest: string }) {
  const aliases = useAliases({ genome_digest: digest, pageSize: 20 });
  const names = aliases.data?.results.map((a) => a.name) ?? [];
  const label = preferredAlias(names) ?? formatDigest(digest);
  return (
    <Link className="rg-link" to={`/genomes/${encodeURIComponent(digest)}?tab=details`}>
      {label}
    </Link>
  );
}
