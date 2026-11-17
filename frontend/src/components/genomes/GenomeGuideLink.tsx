import { ExternalLink } from '../common/ExternalLink';
import { genomeGuide } from '../../services/projectLinks';

/** A link to the docs guide for choosing among genomes like these, if there is one. */
export function GenomeGuideLink({
  terms,
}: {
  terms: readonly (string | null | undefined)[];
}) {
  const guide = genomeGuide(terms);
  if (!guide) return null;
  return (
    <p className="text-sm">
      Not sure which genome to use? Read{' '}
      <ExternalLink href={guide.href}>{guide.title}</ExternalLink>
    </p>
  );
}
