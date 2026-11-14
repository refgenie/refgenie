import { ExternalLink } from '../common/ExternalLink';
import { Badge } from '../common/Badge';
import { accessionUrl } from '../../utils/links';
import { namingLabel, namingTitle } from '../../utils/namingStyle';
import type { Naming } from '../../utils/namingStyle';
import type { GenomeDetailResponse } from '../../types/api';

export interface GenomeHeaderProps {
  /** Undefined while loading: every line shows a placeholder of the same height. */
  genome: GenomeDetailResponse | undefined;
  /** How its sequences are named (`chr1` or `1`), once they load. */
  naming?: Naming | null;
}

function Placeholder({ width }: { width: 'sm' | 'md' | 'lg' }) {
  return (
    <span className={`rg-skeleton rg-skeleton--inline rg-genome-head__ph--${width}`} />
  );
}

/**
 * The key facts under a genome's title: species, accession and description on
 * the left. The digest is in the head's corner; counts and size are on the tabs.
 * Aliases are listed in full on the Details tab.
 */
export function GenomeHeader({ genome, naming }: GenomeHeaderProps) {
  return (
    <div className="rg-genome-head" aria-busy={genome ? undefined : true}>
      <div className="rg-genome-head__facts">
      <p className="rg-genome-head__line">
        {genome ? <Subtitle genome={genome} naming={naming} /> : <Placeholder width="lg" />}
      </p>
      {genome ? (
        genome.description && <Description text={genome.description} />
      ) : (
        <p className="rg-genome-head__line">
          <Placeholder width="lg" />
        </p>
      )}
      </div>
    </div>
  );
}

function Subtitle({
  genome,
  naming,
}: {
  genome: GenomeDetailResponse;
  naming?: Naming | null;
}) {
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
        <ExternalLink key="accession" href={href} className="rg-link rg-link--quiet">
          {accession}
        </ExternalLink>
      ) : (
        <span key="accession" className="rg-muted">
          {accession}
        </span>
      ),
    );
  }
  if (parts.length === 0) {
    parts.push(
      <span key="none" className="rg-muted">
        Species not recorded
      </span>,
    );
  }
  if (naming) {
    parts.push(
      <Badge key="naming" title={namingTitle(naming)}>
        {namingLabel(naming)}
      </Badge>,
    );
  }
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

/** The whole description: no truncation, no "More". */
function Description({ text }: { text: string }) {
  return (
    <div className="rg-genome-head__desc">
      <span className="rg-genome-head__desc-label">Description</span>
      <p>{text}</p>
    </div>
  );
}
