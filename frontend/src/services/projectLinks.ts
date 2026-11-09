/** Links for the refgenie project itself; the same on every instance. */

export const PROJECT_LINKS = {
  docs: 'https://docs.refgenie.org',
  github: 'https://github.com/refgenie/refgenie',
  issues: 'https://github.com/refgenie/refgenie/issues',
  license: 'https://github.com/refgenie/refgenie/blob/master/LICENSE.txt',
  pypi: 'https://pypi.org/project/refgenie/',
  /** Skill file for the CLI and Python API (not this instance's HTTP API). */
  toolSkill: 'https://docs.refgenie.org/SKILL.md',
  useWithAi: 'https://docs.refgenie.org/refgenie/use_with_ai/',
} as const;

export interface Citation {
  readonly authors: string;
  readonly year: string;
  readonly title: string;
  readonly venue: string;
  readonly doi: string;
  readonly note: string;
  /** One-line plain text for the copy button. */
  readonly plain: string;
}

/** From refgenie-docs/docs/refgenie/manuscripts.md. */
export const CITATIONS: readonly Citation[] = [
  {
    authors: 'Stolarczyk M, Reuter VP, Smith JP, Magee NE, Sheffield NC',
    year: '2020',
    title: 'Refgenie: a reference genome resource manager',
    venue: 'GigaScience 9(2)',
    doi: 'https://doi.org/10.1093/gigascience/giz149',
    note: 'The introductory publication. Cite this one if you cite only one.',
    plain:
      'Stolarczyk M, Reuter VP, Smith JP, Magee NE, Sheffield NC. ' +
      'Refgenie: a reference genome resource manager. GigaScience. 2020;9(2). ' +
      'doi:10.1093/gigascience/giz149',
  },
  {
    authors: 'Stolarczyk M, Xue B, Sheffield NC',
    year: '2021',
    title: 'Identity and compatibility of reference genome resources',
    venue: 'NAR Genomics and Bioinformatics 3(2)',
    doi: 'https://doi.org/10.1093/nargab/lqab036',
    note: 'Sequence-derived genome identifiers and provenance tracking, the work behind the digests this interface shows.',
    plain:
      'Stolarczyk M, Xue B, Sheffield NC. Identity and compatibility of ' +
      'reference genome resources. NAR Genomics and Bioinformatics. 2021;3(2). ' +
      'doi:10.1093/nargab/lqab036',
  },
];
