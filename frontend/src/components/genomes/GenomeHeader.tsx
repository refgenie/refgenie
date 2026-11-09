import { DigestChip } from '../common/DigestChip';
import { ExternalLink } from '../common/ExternalLink';
import { accessionUrl } from '../../utils/links';
import { formatBp } from '../../utils/format';
import type { GenomeDetailResponse } from '../../types/api';

export interface GenomeHeaderProps {
  /** Undefined while loading: every line shows a placeholder of the same height. */
  genome: GenomeDetailResponse | undefined;
  /** Null while loading; `false` when the genome has no sequence data. */
  sequenceStats: { count: number; bases: number } | false | null;
  assetCount: number | null;
}

function Placeholder({ width }: { width: 'sm' | 'md' | 'lg' }) {
  return (
    <span className={`rg-skeleton rg-skeleton--inline rg-genome-head__ph--${width}`} />
  );
}

/**
 * The key facts under a genome's title: species, size, digest, description.
 * Aliases are listed in full on the Details tab.
 */
export function GenomeHeader({
  genome,
  sequenceStats,
  assetCount,
}: GenomeHeaderProps) {
  return (
    <div className="rg-genome-head" aria-busy={genome ? undefined : true}>
      <p className="rg-genome-head__line">
        {genome ? <Subtitle genome={genome} /> : <Placeholder width="lg" />}
      </p>
      <p className="rg-genome-head__line rg-genome-head__stats">
        <Stats sequenceStats={sequenceStats} assetCount={assetCount} />
      </p>
      <div className="rg-genome-head__line rg-genome-head__ids">
        {genome ? (
          <DigestChip digest={genome.digest} length={32} />
        ) : (
          <Placeholder width="md" />
        )}
      </div>
      {genome ? (
        genome.description && <Description text={genome.description} />
      ) : (
        <p className="rg-genome-head__line">
          <Placeholder width="lg" />
        </p>
      )}
    </div>
  );
}

function Subtitle({ genome }: { genome: GenomeDetailResponse }) {
  const accession = genome.assembly_accession;
  const href = accessionUrl(
    accession,
    typeof genome.fhr?.accessionID?.url === 'string'
      ? genome.fhr.accessionID.url
      : null,
  );
  const parts = [];
  if (genome.species_name) {
    parts.push(
      <span key="species">
        <em>{genome.species_name}</em>
        {genome.common_name && ` (${genome.common_name})`}
      </span>,
    );
  }
  if (accession) {
    parts.push(
      href ? (
        <ExternalLink key="accession" href={href}>
          {accession}
        </ExternalLink>
      ) : (
        <span key="accession">{accession}</span>
      ),
    );
  }
  if (parts.length === 0) return <span className="rg-muted">Species not recorded</span>;
  return (
    <>
      {parts.map((part, i) => (
        <span key={i}>
          {i > 0 && <span aria-hidden="true"> · </span>}
          {part}
        </span>
      ))}
    </>
  );
}

function Stats({
  sequenceStats,
  assetCount,
}: Pick<GenomeHeaderProps, 'sequenceStats' | 'assetCount'>) {
  const assets =
    assetCount === null ? (
      <Placeholder width="sm" />
    ) : assetCount === 0 ? (
      'no built assets'
    ) : (
      `${assetCount.toLocaleString()} built ${assetCount === 1 ? 'asset' : 'assets'}`
    );
  if (sequenceStats === false) return <>{assets}</>;
  return (
    <>
      {sequenceStats ? (
        `${sequenceStats.count.toLocaleString()} ${sequenceStats.count === 1 ? 'sequence' : 'sequences'} · ${formatBp(sequenceStats.bases)}`
      ) : (
        <Placeholder width="md" />
      )}
      <span aria-hidden="true"> · </span>
      {assets}
    </>
  );
}

/** The whole description: no truncation, no "More". */
function Description({ text }: { text: string }) {
  return <p className="rg-genome-head__desc">{text}</p>;
}
